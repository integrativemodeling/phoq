#!/usr/bin/env python

import unittest
import sys
import shutil
import os
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))

class Tests(unittest.TestCase):
    def test_make_ca_angle(self):
        """Test construction of CA angle restraint"""
        os.chdir(os.path.join(TOPDIR, 'modeling', 'data'))
        with open('file.list', 'w') as fh:
            fh.write(os.path.join(TOPDIR, 'test', '101m.dssp') + '\n')
        p = subprocess.check_call(['./angle.py'])
        # Make sure output was generated
        with open('output.angle.dat') as fh:
            lines = fh.readlines()
        self.assertEqual(len(lines), 58)
        self.assertEqual(lines[1].rstrip('\r\n '), '  ---    86              1')
        os.unlink('output.angle.dat')
        os.unlink('file.list')

    def test_make_ca_dihedral(self):
        """Test construction of CA dihedral restraint"""
        os.chdir(os.path.join(TOPDIR, 'modeling', 'data'))
        with open('file.list', 'w') as fh:
            fh.write(os.path.join(TOPDIR, 'test', '101m.dssp') + '\n')
        p = subprocess.check_call(['./dihedral.py'])
        # Make sure output was generated
        with open('output.dat') as fh:
            lines = fh.readlines()
        self.assertEqual(len(lines), 73)
        self.assertEqual(lines[1].rstrip('\r\n '),
                         '-----    10   100              1')
        os.unlink('output.dat')
        os.unlink('file.list')

    def run_modeling(self, dirname):
        os.chdir(os.path.join(TOPDIR, 'modeling', dirname))
        p = subprocess.check_call(['./HK_model_ENSEMBLE_REM.py', '--test'])
        # todo: make sure the outputs are actually reasonable
        os.unlink('log0')
        os.unlink('traj0.rmf')

    def test_1_state(self):
        """Test 1-state modeling"""
        self.run_modeling('2Y20_MC_1_FULL')

    def test_2_state(self):
        """Test 2-state modeling"""
        self.run_modeling('2Y20_MC_2_FULL')

    def test_2_state_90_1(self):
        """Test 2-state modeling with 90% data (first set)"""
        self.run_modeling('2Y20_MC_2_FULL_90_PERCENT_DATA_1')

    def test_2_state_90_2(self):
        """Test 2-state modeling with 90% data (second set)"""
        self.run_modeling('2Y20_MC_2_FULL_90_PERCENT_DATA_2')

    def test_2_state_95_1(self):
        """Test 2-state modeling with 95% data (first set)"""
        self.run_modeling('2Y20_MC_2_FULL_95_PERCENT_DATA_1')

    def test_2_state_95_2(self):
        """Test 2-state modeling with 95% data (second set)"""
        self.run_modeling('2Y20_MC_2_FULL_95_PERCENT_DATA_2')

    def test_3_state(self):
        """Test 3-state modeling"""
        self.run_modeling('2Y20_MC_3_FULL')

if __name__ == '__main__':
    unittest.main()
