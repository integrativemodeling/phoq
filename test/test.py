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

if __name__ == '__main__':
    unittest.main()
