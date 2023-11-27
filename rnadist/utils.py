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

def list_internal_leafsets(node):

    # if is leaf
    if len(node.children)==0:
        return []

    # local internal_leaf set
    internal_leafset = []
    for c in node.children:
        if len(c.children)==0:
            internal_leafset.append(c.label)

    internal_leafsets = [internal_leafset]

    # adding leafsets of the children
    for c in node.children:
        internal_leafsets += list_internal_leafsets(c)

    return internal_leafsets

def fitch_with_trait_vector(phylo_T, vector_traits):
    """
    Input:
        - phylo_T: the phylogenetic tree. Each node has a name/label. 
        The labels of the leaves are integers.
        - vector_traits: a dictionnary associating each leaf label with
        a vector of traits (numpy vector of zeros and ones). This dictionary
        will be augmented with new assignments and returned at the end.
    """

    # retrieving N, the total number of traits
    N = vector_traits.N

    # B: the "B" variable of the Fitch algorithm
    B = {}
    for c in range(N):
        B[c] = {}

    # auxiliary recursive function to fill B dictionary
    def fillB(c, node):
        if node.isleaf:
            B[c][node.label] = set([vector_traits[node.label][c]])
            return B[c][node.label]
        else:
            v, w = node.children
            Bv = fillB(c, v)
            Bw = fillB(c, w)
            if len(Bv.intersection(Bw)) == 0:
                B[c][node.label] = Bv.union(Bw)
            else:
                B[c][node.label] = Bv.intersection(Bw)

            return B[c][node.label]

    # calling it for each feature
    for c in range(N):
        fillB(c, phylo_T)

    # F; the final assignment of Fitch
    F = {}
    for c in range(N):
        F[c] = {}

    # auxiliary recursive function to fill F
    def fillF(c, node):

        for child in node.children:
            if F[c][node.label] in B[c][child.label]:
                F[c][child.label] = F[c][node.label]
            else:
                F[c][child.label] = B[c][child.label][0]
            fillF(c, child)

    # for each trait, fill F for root and call on root
    for c in range(N):
        if 0 in B[c][phylo_T.label]:
            F[c][phylo_T.label] = 0
        else:
            F[c][phylo_T.label] = 1
        fillF(c, phylo_T)

    # final annotation of the tree (pure interface issue, nothing is computed)
    def build_vector_traits(node):
        if not node.isleaf:
            vector = np.zeros(N)
            for c in range(N):
                vector[c] = F[c][node.label]
            vector_traits[node.label] = vector
            for child in node.children:
                build_vector_traits(child)

    # call for annotation
    build_vector_traits(phylo_T)

    return vector_traits

class NewickNode:
    def __init__(self, name):
        self.descendants = []
        self.name = name

def parse_newick(newick_string):
    if newick_string.endswith(';'):
        root = NewickNode(root)


