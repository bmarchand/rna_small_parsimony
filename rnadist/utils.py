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
        sub_clades.append(frozenset(sorted(clade)))

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

def gapless_aligned_structs_to_clade_traits():
    """
    
    """
    pass

class PhyloNode:

    def __init__(self, label, isleaf, annotation=None):
        self.isleaf = isleaf
        self.label = label
        self.annotation = annotation
        self.children = []

    def add_children(self, node1, node2):
        self.children = [node1, node2]


def tree_tab_file_parse(tree_tab_lines):

    if len(tree_tab_lines) > 1:
    # not a leaf
        first_line = tree_tab_lines[0]
        root_node = PhyloNode(first_line.split(' ')[0].rstrip('\n'), False)
        if first_line.find(':') >=0:
            root_node.annotation = first_line.split(' ')[-1].rstrip('\n')

        # children
        child1_lines = [tree_tab_lines[1]]
        index = 2
        while tree_tab_lines[index].startswith('\t\t'):
            child1_lines.append(tree_tab_lines[index])
            index += 1
        child2_lines = [tree_tab_lines[index]]
        index += 1
        while index < len(tree_tab_lines) and tree_tab_lines[index].startswith('\t\t'):
            child2_lines.append(tree_tab_lines[index])
            index += 1

        child1_lines = [l[1:] for l in child1_lines]
        child2_lines = [l[1:] for l in child2_lines]

        child1 = tree_tab_file_parse(child1_lines)
        child2 = tree_tab_file_parse(child2_lines)
        
        root_node.add_children(child1, child2)

    else:
        first_line = tree_tab_lines[0]
        root_node = PhyloNode(first_line.split(' ')[0].rstrip('\n'), True)
        if first_line.find(':') >=0:
            root_node.annotation = first_line.split(' ')[-1].rstrip('\n')

    return root_node

def list_leaves(node):
    if node.isleaf:
        return [node.label]
    else:
        L = []
        for child in node.children:
            L += list_leaves(child)
        return L

def depthk_bps(structure, depth=0):

    L = []
    stack = []

    for k, c in enumerate(structure):
        if c=='(':
            stack.append(k)
        if c==')':
            l = stack.pop()
            if len(stack)==depth:
                L.append((l,k))

    return L

def build_tree(structure):

    L0 = depthk_bps(structure)
    assert(len(L0)==1)

    root = Node()
    i,j = L0[0]
    Li = Node(label=i)
    Lj = Node(label=j)
    root.add_child(Li)
    root.add_child(Lj)

    # rest of the leaf set
    IL = list(range(i+1,j,1))

    L1 = depthk_bps(structure, depth=1)

    for k,l in L1:
        assert(structure[k]=='(')
        assert(structure[l]==')')
        IL = [m for m in IL if m < k or m > l]
        # creating some leaves
        for m in IL:
            if m < k:
                Lm = Node(label=m)
                root.add_child(Lm)

        sub_structure = list(structure)
        for index in range(len(sub_structure)):
            if index < k or index > l:
                sub_structure[index] = '.'
        child_kl = build_tree(''.join(sub_structure))

        root.add_child(child_kl)

    # creating the rest of the leaves
    if len(L1) > 0:
        M = max([max(k,l) for k,l in L1])
        for m in IL:
            if m > M:
                Lm = Node(label=m)
                root.add_child(Lm)
    else:
        for m in IL:
            Lm = Node(label=m)
            root.add_child(Lm)


    return root
    
def check_well_parenthesized(s):
    stack = []
    for k,c in enumerate(s):
        if c=='(':
            stack.append(k)
        if c==')':
            stack.pop()
    assert(len(stack)==0)

def aligned_gapless_structures_to_trait_vectors(str_dict):

    lengths = set([len(val) for _, val in str_dict.items()])
    assert(len(lengths)==1)
    length = lengths.pop()

    tree_dict = {}

    trait_vector = VectorTraits(length)

    clade_set = set([])            
    for key, value in str_dict.items():
        value = '('+value+')' # virtual overlapping arc 
        check_well_parenthesized(value)
        tree_dict[key] = build_tree(value)
        for c in list_clades(tree_dict[key]):
            clade_set.add(frozenset(c))
        
    clade_set = list(clade_set) # order is random here

    for key, value in tree_dict.items():
        vector = []
        for k, C in enumerate(clade_set):
            if C in list_clades(value):
                vector.append(1)
            else:
                vector.append(0)

        trait_vector.traits[key] = vector

    return trait_vector, clade_set

def trait_vector_to_str(trait_vector, clade_set, N):

    structure = list('.'*(N+2))

    for k, val in enumerate(trait_vector):
        if val==1:
            i = min(clade_set[k])
            j = max(clade_set[k])
            structure[i] = '('
            structure[j] = ')'

    return ''.join(structure[1:-1])


class VectorTraits:

    def __init__(self, N):
        self.N = N
        self.traits = {}

    def add_trait(self, leaf, vector):
        self.traits[leaf] = vector

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
            B[c][node.label] = set([vector_traits.traits[node.label][c]])
            return B[c][node.label]
        else:
            v, w = tuple(node.children)
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
                F[c][child.label] = list(B[c][child.label])[0]
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
            vector = [0 for _ in range(N)]
            for c in range(N):
                vector[c] = F[c][node.label]
            vector_traits.traits[node.label] = vector
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


