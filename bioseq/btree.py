

class BTree:
    def __init__(self, id=None, left=None, right=None):
        self.id = id
        self.left = left
        self.right = right
        self.height = 0

    def merge(n1, n2, distance):
        n = BTree(left=n1, right=n2)
        d = distance / 2
        n1.height = d
        n2.height = d
        return n

    def to_clade(self):
        if self.left: return '<clade branch_length="%.3f">' % self.height + self.left.to_clade() + self.right.to_clade() + '</clade>'
        return '<clade branch_length="%.3f"><name>%s</name></clade>' % (self.height, self.id)

    def to_xml(self):
        return '<?xml version="1.0" encoding="UTF-8"?><phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd" xmlns="http://www.phyloxml.org"><phylogeny rooted="true"><name>ABI Tree</name><clade>' + (self.left.to_clade() if self.left else "") + (self.right.to_clade() if self.right else "") + "</clade></phylogeny></phyloxml>"
