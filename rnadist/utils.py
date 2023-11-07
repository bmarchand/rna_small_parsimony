class Node:

    def __init__(self, label=None):
        self.children = []
        self.label = label

    def add_child(self, node):
        self.children.append(node)


def RFdistance(root1, root2):
    L1 = list_clades(root1)
    L2 = list_clades(root2)

    unique1 = [c for c in L1 if c not in L2]
    unique2 = [c for c in L2 if c not in L1]

    return len(unique1)+len(unique2)

def list_clades(root):
    """
    Input: 
        - an instance of the Node class, typically the root of a tree
    Output:
        - the list of all clades of the tree. Including the trivial one.
    """
    
    def clade_and_subclades(node):
        """
        Auxiliary recursive function

        Input:
            - instance of Node
        Output:
            - a tuple (L,C) with:
                * C the clade associated to node
                * L a list of all clades below node, including C
        """
        if len(node.children)==0:
            return [[node.label]], [node.label]

        sub_clades = []
        clade = []
        for c in node.children:
            sub, clade_c = clade_and_subclades(c)
            clade += clade_c
            sub_clades += sub
        sub_clades.append(sorted(clade))

        return sub_clades, clade

    # main call
    L, _ = clade_and_subclades(root)
    return L
