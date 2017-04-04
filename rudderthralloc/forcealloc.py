""" forcealloc.py - Map commanded thrust to generalized, cartesian forces """
# BSD 2-Clause License
#
# Copyright (c) 2001-2017, Karl-Petter Lindegaard
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np
from math import sin, cos
from cs2data import T1LX, T1LY, T2LX, T2LY, T3LX, T3LY


class ForceAllocation:
    """
    ForceAllocation - maps from tau_c to generalized forces
    """
    def __init__(self, theta1, theta2, c1, c2):
        """
        :param theta1: Port main propeller/rudder postive force angle span [rad]
        :param theta2: Starboard main propeller/rudder postive force angle span [rad]
        :param c1: Port main propeller positive thrust bias [N]
        :param c2: Starboard main properller positive thrust bias [N]
        """
        self.theta1 = theta1
        self.theta2 = theta2
        self.c1 = c1
        self.c2 = c2

        self.Q1 = np.eye(4)
        self.Q2 = np.eye(4)

        # Full allocation matrix
        self.A = np.array([
            [1, 0, 1, 0, 0],
            [0, 1, 0, 1, 1],
            [-T1LY, T1LX, -T2LY, T2LX, T3LX]
        ])

        # Filters. f1 = Rudder 2 inactive, f2 = Rudder 1 inactive
        self.f1 = np.array([True, True, True, False, True])
        self.f2 = np.array([True, False, True, True, True])

        # Configure A1, n1 and A1dagger etc.
        self.A1 = self.A[:,self.f1]

        # Null-vector for A1
        self.n1 = np.zeros(4)
        self.n1[0] = (T1LX - T3LX) / (T1LY - T2LY)
        self.n1[1] = 1.0
        self.n1[2] = -self.n1[0]
        self.n1[3] = -1.0

        # A1_dagger = A1'*inv(A1*A1')
        self.A1_dagger = self.A1.T.dot(np.linalg.inv(self.A1.dot(self.A1.T)))

        # Configure A2, n2 and A2dagger etc.
        self.A2 = self.A[:,self.f2]

        # Null-vector for A2
        self.n2 = np.zeros(4)
        self.n2[0] = (T3LX - T2LX) / (T1LY - T2LY)
        self.n2[1] = -self.n2[0]
        self.n2[2] = -1.0
        self.n2[3] = 1.0

        # A2_dagger = A2'*inv(A2*A2')
        self.A2_dagger = self.A2.T.dot(np.linalg.inv(self.A2.dot(self.A2.T)))

    def nullsub1(self, tauc, Adagger, n, theta, c):
        # type: (np.array, np.array, np.array, float, float) -> np.array

        # Step 0: Prepare the a-vector (sector boundary)
        a1 = cos(theta)
        a2 = -sin(theta)

        # Step 1: Find optimal solution based on pseudo-inverse
        u0 = Adagger.dot(tauc)

        # Step 2: Extract prop/rudder and translate to "the other" ref. frame
        u0m1 = u0[0] - c
        u0m2 = u0[1]

        # Step 3: Sector check
        nn1 = n[1]
        nn2 = -n[0]
        dp = nn1*u0m1 + nn2*u0m2

        insector = False
        if dp <= 0.0:
            # Traverse in x-asxis (fx,0)
            b1 = 0.0
            b2 = 1.0
        else:
            # Are we in sector "1"
            if u0m2 >= 0.0:
                b1 = 0.0
                b2 = 1.0

            # Or perhaps we are already within the valid sector
            elif u0m1*a2 < u0m2*a1:
                insector = True

            # Otherwise, traverse along the nullvector until sector limit "a"
            else:
                b1 = a2
                b2 = -a1

        # Step 4: Find lambda, the distance to traverse
        gamma = 0.0
        if not insector:
            gamma = -(u0m1*b1 + u0m2*b2) / (n[0]*b1 + n[1]*b2)

        # Step 5: Adjust solution
        u = u0 + gamma*n
        return u

    def nullsub2(self, tauc, Adagger, n, theta, c):
        # type: (np.array, np.array, np.array, float, float) -> np.array

        # Step 0: Prepare the a-vector (sector boundary)
        a1 = cos(theta)
        a2 = sin(theta)

        # Step 1: Find optimal solution based on pseudo-inverse
        u0 = Adagger.dot(tauc)

        # Step 2: Extract prop/rudder and translate to "the other" ref. frame
        u0m1 = u0[1] - c
        u0m2 = u0[2]

        # Step 3: Sector check
        nn1 = n[2]
        nn2 = -n[1]
        dp = nn1 * u0m1 + nn2 * u0m2

        insector = False
        if dp >= 0.0:
            # Traverse in x-asxis (fx,0)
            b1 = 0.0
            b2 = 1.0
        else:
            # Are we in sector "1"
            if u0m2 <= 0.0:
                b1 = 0.0
                b2 = 1.0

            # Or perhaps we are already within the valid sector
            elif u0m1 * a2 > u0m2 * a1:
                insector = True

            # Otherwise, traverse along the nullvector until sector limit "a"
            else:
                b1 = a2
                b2 = -a1

        # Step 4: Find lambda, the distance to traverse
        gamma = 0.0
        if not insector:
            gamma = -(u0m1 * b1 + u0m2 * b2) / (n[1] * b1 + n[2] * b2)

        # Step 5: Adjust solution
        u = u0 + gamma * n
        return u

    def allocate(self, tau):
        """
        Map 3-DOF commanded thrust to generalized forces. First two elements are surge and sway
        for thruster 1 (port main prop+rudder), next two for starboard main prop+rudder, fifth
        element is the bow thruster's sway force.
        :param tau: Commanded thrust vector (surge, sway, yaw)
        :return: Generalized forces
        """
        # type: (np.array) -> np.array

        # Call subroutines for each rudder
        x1 = self.nullsub1(tau, self.A1_dagger, self.n1, self.theta1, self.c1)
        x2 = self.nullsub2(tau, self.A2_dagger, self.n2, self.theta2, self.c2)

        # Compare results and pick the best solution J = x'*Q*x
        j1 = x1.dot(self.Q1.dot(x1))
        j2 = x2.dot(self.Q2.dot(x2))
        u = np.zeros(5)
        if j1 <= j2:
            # Use: u = [x1(0) x1(1) x1(2) 0 x1(3)];
            u[self.f1] = x1
        else:
            # u = [x2(0) 0 x2(1) x2(2) x2(3)];
            u[self.f2] = x2
        return u
