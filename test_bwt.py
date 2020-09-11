# -*- coding: utf-8 -*-
""" Tests for the paleni-bwt module """

import unittest
import bwt as bwt

class TestBWT(unittest.TestCase):

    def test_transform(self):
        ## check transform
        self.assertEqual(bwt.bwTransform("ciao"), 'oi$ca')
        self.assertEqual(bwt.bwTransform("itopinon"), 'np$ointoi')

    def test_inverse(self):
        ## check inverse
        self.assertEqual(bwt.bwInverse("oi$ca"), 'ciao')
        self.assertEqual(bwt.bwInverse("np$ointoi"), 'itopinon')

    def test_type_error(self):
        s = 56
        ## check that fails when the input is not a string
        with self.assertRaises(TypeError):
            bwt.bwTransform(s)
        with self.assertRaises(TypeError):
            bwt.bwInverse(s)
    
    def test_value_error(self):
        ## check that fails if given a string with wrong format
        s = "ciao$"
        s2 = "ciao"
        s3 = "ci$ao$"
        
        with self.assertRaises(bwt.TerminatorError) as e1:
            bwt.bwTransform(s)
        self.assertEqual(e1.exception.error_code, 1)
        with self.assertRaises(bwt.TerminatorError) as e2:
            bwt.bwInverse(s2)
        self.assertEqual(e2.exception.error_code, 2)
        with self.assertRaises(bwt.TerminatorError) as e3:
            bwt.bwInverse(s3)
        self.assertEqual(e3.exception.error_code, 3)
    
    

    def test_alternative_terminators(self):
        s ="oi%oac"
        with self.assertRaises(bwt.TerminatorError) as e1:
            bwt.bwInverse(s)
        self.assertEqual(e1.exception.error_code, 2)
        with self.assertRaises(bwt.TerminatorError) as e2:
            bwt.bwTransform(s, 2)
        self.assertEqual(e2.exception.error_code, 4)
        with self.assertRaises(bwt.TerminatorError) as e3:
            bwt.bwTransform(s, "%$")
        self.assertEqual(e3.exception.error_code, 5)
        with self.assertRaises(bwt.TerminatorError) as e4:
            bwt.bwTransform(s, "+")
        self.assertEqual(e4.exception.error_code, 6)
        
        

if __name__ == '__main__':
    unittest.main()