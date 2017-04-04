""" test_forcealloc.py - unit tests for forcealloc.py """
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

import unittest
import numpy as np
from numpy.linalg import norm
from forcealloc import ForceAllocation


class TestForceAllocation(unittest.TestCase):

    def test_nullvectors(self):
        fa = ForceAllocation(0, 0, 0, 0)
        test_1 = fa.A1.dot(fa.n1)
        test_2 = fa.A2.dot(fa.n2)
        self.assertAlmostEqual(norm(test_1), 0, 8, msg="A1*n1 != 0")
        self.assertAlmostEqual(norm(test_2), 0, 8, msg="A2*n2 != 0")

    def test_inverses(self):
        fa = ForceAllocation(0, 0, 0, 0)
        test_1 = fa.A1.dot(fa.A1_dagger) - np.eye(3)
        test_2 = fa.A2.dot(fa.A2_dagger) - np.eye(3)
        self.assertAlmostEqual(norm(test_1), 0.0, 8, msg="A1*A1_dagger != I")
        self.assertAlmostEqual(norm(test_2), 0.0, 8, msg="A2*A2_dagger != I")

    def test_surge(self):
        fa = ForceAllocation(0, 0, 0, 0)
        u = fa.allocate([10, 0, 0])
        self.assertAlmostEqual(u[0], 5.0, 10)
        self.assertAlmostEqual(u[2], 5.0, 10)

    def test_yaw_no_rudder(self):
        fa = ForceAllocation(0, 0, 0, 0)
        u = fa.allocate([0, 0, 10])
        self.assertAlmostEqual(u[0], -u[2], 10)
        self.assertAlmostEqual(u[1], 0.0, 10)
        self.assertAlmostEqual(u[3], 0.0, 10)
        self.assertAlmostEqual(u[4], 0.0, 10)
        tau = fa.A.dot(u)
        self.assertAlmostEqual(tau[0], 0.0, 10)
        self.assertAlmostEqual(tau[1], 0.0, 10)
        self.assertAlmostEqual(tau[2], 10.0, 10)
