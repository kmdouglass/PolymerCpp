# © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE,
# Switzerland, Laboratory of Experimental Biophysics, 2017
# See the LICENSE.txt file for more details.

import unittest
from PolymerCpp.helpers import getCppWLC
    
class ChainGenerationTest(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_getCppWLC(self):
        """Basic unittest for method getCppWLC.

        """
        pathLength = 50
        chain      = getCppWLC(pathLength=pathLength)

        self.assertEqual(chain.shape, (51, 3))
    
if __name__ == '__main__':
    unittest.main()
