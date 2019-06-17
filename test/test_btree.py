import unittest
from bioseq import *

class TestBTree(unittest.TestCase):
    def test_constructor(self):
        b = BTree()
        self.assertIsNone(b.id)
        self.assertIsNone(b.left)
        self.assertIsNone(b.right)
        self.assertEqual(b.height, 0)

        b1 = BTree("1", b, b)
        self.assertEqual(b1.id, "1")
        self.assertEqual(b1.left, b)
        self.assertEqual(b1.right, b)

    def test_merge(self):
        b = BTree.merge(BTree("1"), BTree("2"), 10)
        self.assertEqual(b.height, 0)
        self.assertEqual(b.left.id, "1")
        self.assertEqual(b.left.height, 5.0)
        self.assertEqual(b.right.height, 5.0)

    def test_to_clade(self):
        b = BTree.merge(BTree("1"), BTree("2"), 10)
        self.assertEqual(b.to_clade(), '<clade branch_length="0.000"><clade branch_length="5.000"><name>1</name></clade><clade branch_length="5.000"><name>2</name></clade></clade>')

    def test_to_xml(self):
        b = BTree.merge(BTree("1"), BTree("2"), 10)
        self.assertEqual(b.to_xml(), '<?xml version="1.0" encoding="UTF-8"?><phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd" xmlns="http://www.phyloxml.org"><phylogeny rooted="true"><name>ABI Tree</name><clade><clade branch_length="5.000"><name>1</name></clade><clade branch_length="5.000"><name>2</name></clade></clade></phylogeny></phyloxml>')
