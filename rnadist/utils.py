import numpy as np

class Node:

    def __init__(self, label=None):
        self.children = []
        self.label = label

    def add_child(self, node):
        self.children.append(node)

class PhyloNode:

    def __init__(self, label, isleaf, annotation=None):
        self.isleaf = isleaf
        self.label = label
        self.annotation = annotation
        self.children = []

    def add_children(self, node1, node2):
        self.children = [node1, node2]

def pairing_info(ss_cons):
    '''
    Getting the pairing information of a structure
    
    Input:
        - structure as string in dot bracket notation

    Output:
        - d: dictionary which associates each paired
        position to its partner
        - bps: the list of base pairs
    '''
    
    parenthesis_systems = ['()','[]','{}','<>']

    stacks = {}
    for oc in parenthesis_systems:
        stacks[oc] = []

    d = {}
    bps = []

    for k, c in enumerate(ss_cons):
        for oc in parenthesis_systems:
            if c==oc[0]:
                stacks[oc].append(k)
            if c==oc[1]:
                l = stacks[oc].pop()
                d[k] = l
                d[l] = k 
                bps.append((k,l))

    return d, bps


def RFdistance(root1, root2):
    L1 = list_clades(root1)
    L2 = list_clades(root2)

    unique1 = [c for c in L1 if c not in L2]
    unique2 = [c for c in L2 if c not in L1]

    return len(unique1)+len(unique2)



def C2_ILMedian(input_structures):

    # adding virtual overarching base pair
    input_structures = ['('+s+')' for s in input_structures]

    # number of positions
    N = len(input_structures[0])

    # base pairs
    bps = []
    for s in input_structures:
        assert(len(s)==N)
        _, l = pairing_info(s)
        bps += l

    # input leaf sets
    input_leaf_sets = []
    input_leaf_set_dict = {}
    for s in input_structures:
        tree = build_tree(s)
        ILs = list_internal_leafsets(tree)
        input_leaf_sets += ILs
        input_leaf_set_dict[s] = ILs 

    # auxiliary recursive function 1: optimal score
    print(input_leaf_sets)
    print(N)

    def optimal_score(i,j):
        print("call",i,j)
        if (i,j) in C.keys():
            print('-> already computed')
            return C[(i,j)]

        C[(i,j)] = np.inf
        for I in input_leaf_sets:
            if min(I)==i and max(I)==j:
                print('candidate I', I)
                score_with_I = len([key for key, val in input_leaf_set_dict.items() if I in val])
                holes = IL_holes(I)
                for i_h, j_h in holes:
                    score_with_I += optimal_score(i_h+1,j_h-1)

                if score_with_I < C[(i,j)]:
                    C[(i,j)] = score_with_I

        return C[(i,j)]

    # auxiliary recursive function 2: backtrace
    def backtrace(i,j):
        for I in input_leaf_sets:
            if min(I)==i and max(I)==j:
                print(i,j,':','candidate leaft set', I)
                score_with_I = len([key for key, val in input_leaf_set_dict.items() if I in val])
                holes = IL_holes(I)
                print('score with & holes',score_with_I, holes)
                for i_h, j_h in holes:
                    score_with_I += C[(i_h+1,j_h-1)]

                if score_with_I==C[(i,j)]:
                    result = [I]
                    for i_h, j_h in IL_holes(I):
                        result += backtrace(i_h+1, j_h-1)
                    return result

    
    C = {} # DP table

    OPT = optimal_score(0, N-1)
    print(C)
        
    ILs = backtrace(0,N-1)

    median = list('.'*N)

    for I in ILs:
        print('I from final ILs list', I)
        median[min(I)] = '('
        median[max(I)] = ')'

    return ''.join(median[1:-1])

def IL_holes(I):
    I = sorted(I)
    holes = []
    for k in range(len(I)-1):
        if I[k+1]>I[k]+1:
            holes.append((I[k],I[k+1]))
    return holes

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
            return [frozenset([node.label])], [node.label]

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


def tree_tab_file_parse(tree_tab_lines):

    if len(tree_tab_lines) > 1:
    # not a leaf
        first_line = tree_tab_lines[0]
        root_node = PhyloNode(first_line.split(' ')[0].rstrip('\n'), False)
        if first_line.find(':') >=0:
            root_node.annotation = first_line.split(' ')[-1].rstrip('\n')

        # children
        children_lines = {}

        for k, line in enumerate(tree_tab_lines[1:]):
            if not line.startswith('\t\t'):
                children_lines[line] = []
                current_line = line

            if line.startswith('\t\t'):
                children_lines[current_line].append(line)

        children_lines = [[key]+value for key, value in sorted(children_lines.items(), key=lambda u: u[0])]
        children_lines = [[l[1:] for l in child] for child in children_lines]

        for child in children_lines:
            root_node.children.append(tree_tab_file_parse(child))

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


def median_based_heuristic(phylo_T, 
                           str_dict, 
                           rounds=10, 
                           until_converge=False, 
                           median_function=C2_ILMedian):
    """
        - phylo_T: phylogenetic tree.
        - str_dict: dictionary label -> str. only on leaves at first, to complete
        and return
    """
    
    random_key = np.random.choice(str_dict.keys())
    random_input_str = str_dict[random_key]

    queue = [phylo_T]
    while len(queue) > 0:
        node = queue.pop()
        str_dict[node.label] = random_input_str

    cnt = 0
    keep_going = True
    while keep_going:
        keep_going = False
        cnt += 1

        queue = [(phylo_T, c) for c in phylo_T.children]

        smth_changed = False
        while len(queue) > 0 and not smth_changed:
            parent, child = queue.pop()
            if len(child.children) > 0:
                # current sum over edges around node
                current_cost = IL_distance(str_dict[parent.label], 
                                           str_dict[child.label])
                for grand_child in child.children:
                    current_cost += IL_distance(str_dict[child.label],
                                                str_dict[grand_child.label])

                # median computation
                neighbor_list = [parent] + child.children
                str_list = [str_dict[n.label] for n in neighbor_list]
                median = median_function(str_list)

                # new sum over edges around node
                new_cost = IL_distance(str_dict[parent.label], 
                                       median)
                for grand_child in child.children:
                    new_cost += IL_distance(str_dict[child.label],
                                                str_dict[grand_child.label])
                
                # if new cost better, change label of node to median
                if new_cost < current_cost:
                    str_dict[child.label] = median
                    smth_changed = True

        # do we keep going ? count condition, and whether smth changed
        if cnt < rounds and smth_changed:
            keep_going = True

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
            num_1 = 0
            num_0 = 0
            for v in node.children:
                Bc = fillB(c, v)
                if Bc==set([1]):
                    num_1 += 1
                if Bc==set([0]):
                    num_0 += 1
            if num_0 > num_1:
                B[c][node.label] = set([0])
            elif num_0 < num_1:
                B[c][node.label] = set([1])
            elif num_0==num_1:
                B[c][node.label] = set([0,1])

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


