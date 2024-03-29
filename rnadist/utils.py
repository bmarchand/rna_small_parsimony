import numpy as np
import re
from rnadist.indset_intgraph import get_max_interval_indset

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


def RFdistance(str1, str2):

    root1 = build_tree('('+str1+')')
    root2 = build_tree('('+str2+')')
    
    L1 = list_clades(root1)
    L2 = list_clades(root2)

    unique1 = [c for c in L1 if c not in L2]
    unique2 = [c for c in L2 if c not in L1]

    return len(unique1)+len(unique2)

def C2_ILMedian(input_structures):
    return median_structure(input_structures)

def C2_RFMedian(input_structures):
    return median_structure(input_structures, metric='RF')

def unconstrained_ILMedian(input_structures):
    print(input_structures)
    return median_structure(input_structures, input_leafsets_only=False)

def unconstrained_RFMedian(input_structures):
    
    count_bps = {}

    for structure in input_structures:
        d, bps = pairing_info(structure)

        for bp in bps:
            try:
                count_bps[bp] += 1
            except KeyError:
                count_bps[bp] = 0

    N = len(structure[0])
    median = list('.' * N)

    for bp, cnt in count_bps.items():
        if cnt > len(input_structures)/2:
            median[bp[0]] = '('
            median[bp[1]] = ')'

    return median

def filter_dot_only(input_structures):

    dot_positions = []

    N = len(input_structures[0])

    for k in range(N):
        all_points = True
        for structure in input_structures:
            if structure[k]!='.':
                all_points = False
                break

        if all_points:
            dot_positions.append(k)

    new_structures = []
    for structure in input_structures:
        l = list(structure)
        for k in dot_positions:
            l[k]=''

        new_structures.append(''.join(l))

    return dot_positions, new_structures

def reinsert_dots(median, dot_positions):

    l = list(median)

    for k in sorted(dot_positions, key=lambda x: x):
        l.insert(k,'.')

    return ''.join(l)

def median_structure(input_structures, metric='IL', input_leafsets_only=True):
    """
        General function for computing the median of input RNA structures.

        Input:
            - input_structures: list of input RNA structures in dot-bracket notation.
            - metric: (default value=IL) The distance with respect to which the median
            is computed. Robinson-Foulds if metric=RF, internal leaf-set distance
            if metric=IL
            - input_leafsets (default value=True). Whether or not the search space should
            be restricted to input leaf-sets only, or not.
        
        Output:
            - the median structure in dot-bracket notation.
    """

    dot_positions, input_structures = filter_dot_only(input_structures)
    print('reduced input_structures', input_structures)

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
    for i, s in enumerate(input_structures):
        tree = build_tree(s)
        ILs = list_internal_leafsets(tree)
        input_leaf_sets += ILs
        input_leaf_set_dict[str(i)+'-'+s] = ILs 

    # input clades
    input_clades_dict = {}
    for i, s in enumerate(input_structures):
        tree = build_tree(s)
        clades = list_clades(tree)
        input_clades_dict[str(i)+'-'+s] = clades



    # auxiliary recursive function 1: optimal score

    def optimal_score(i,j):
        if (i,j) in C.keys():
            return C[(i,j)]

        C[(i,j)] = np.inf

        # looping over all input leaf sets
        for I in input_leaf_sets:
            if min(I)==i and max(I)<=j:
                if metric=='IL':
                    score_with_I = len([key for key, val in input_leaf_set_dict.items() if I not in val])
                    score_with_I -= len([key for key, val in input_leaf_set_dict.items() if I in val])
                if metric=='RF':
                    clade_I = frozenset(list(range(min(I),max(I)+1,1)))
                    score_with_I = len([key for key, val in input_clades_dict.items() if I not in val])
                    score_with_I -= len([key for key, val in input_clades_dict.items() if I in val])
                holes = IL_holes(I,i,j)
                for i_h, j_h in holes:
                    score_with_I += optimal_score(i_h,j_h)

                if score_with_I < C[(i,j)]:
                    C[(i,j)] = score_with_I

        if not input_leafsets_only:
            # maximum weighted independent set.
            sol = max_weigthed_is_result(i,j)

            weight = sum( i[2] for i in sol )

            C[(i,j)] = min(C[(i,j)], len(input_structures)-weight)

        return C[(i,j)]

    # auxiliary recursive function 2: backtrace
    def backtrace(i,j):

        # same loop over input leaf-sets
        for I in input_leaf_sets:
            if min(I)==i and max(I)<=j:
                if metric=='IL':
                    score_with_I = len([key for key, val in input_leaf_set_dict.items() if I not in val])
                    score_with_I -= len([key for key, val in input_leaf_set_dict.items() if I in val])
                if metric=='RF':
                    score_with_I = len([key for key, val in input_clades_dict.items() if I not in val])
                    score_with_I -= len([key for key, val in input_clades_dict.items() if I in val])

                holes = IL_holes(I,i,j)
                for i_h, j_h in holes:
                    score_with_I += C[(i_h,j_h)]

                # if the score with I chosen is optimal: I part of the solution.
                if score_with_I==C[(i,j)]:
                    result = [I]
                    for i_h, j_h in IL_holes(I,i,j):
                        result += backtrace(i_h,j_h)
                    return result

        if not input_leafsets_only:
            # maximum weighted independent set.
            sol = max_weigthed_is_result(i,j)

            indset_intervals = [(i[0],i[1]) for i in sol]
            weight = sum( i[2] for i in sol )

            if C[(i,j)]==len(input_structures)-weight:
                I = list(range(i,j+1,1))
                for u,v in indset_intervals:
                    for x in range(u,v+1,1):
                        I.remove(x)

                result = [I]
                for i_h, j_h in IL_holes(I,i,j):
                    result += backtrace(i_h, j_h)
                return result
    
    # auxiliary function 3: maximum weighted independent set on (i,j)
    MWIS = {}
    def max_weigthed_is_result(i,j):
        if (i,j) in MWIS.keys():
            return MWIS[(i,j)]
            
        intervals = []
        for u in range(i+1,j,1):
            for v in range(u+1, j,1):
                intervals.append((u,v,-optimal_score(u,v)))

        sol = get_max_interval_indset(intervals)

        MWIS[(i,j)] = sol
        return MWIS[(i,j)]
    
    C = {} # DP table

    OPT = optimal_score(0, N-1)
        
    ILs = backtrace(0,N-1)

    covered = {}
    for i in range(N):
        covered[i]=False
    for I in ILs:
        for i in I:
            assert(not covered[i])
            covered[i] = True
    for i in range(N):
        assert(covered[i])

    median = list('.'*N)

    for I in ILs:
        median[min(I)] = '('
        median[max(I)] = ')'
        

    check_well_parenthesized(''.join(median[1:-1]))

    median = ''.join(median[1:-1])

    return reinsert_dots(median, dot_positions) 

def IL_holes(I,i,j):
    '''
    Positions of [i,j] not covered by I, regrouped in intervals
    '''
    I = sorted(I)
    holes = []
    for k in range(len(I)-1):
        if I[k+1]>I[k]+1:
            holes.append((I[k]+1,I[k+1]-1))

    if i not in I:
        holes = [(i,I[0]-1)]+holes
    if j not in I:
        holes += [(I[-1]+1,j)]
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
    L = [c for c in L if len(c)>1]
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

def IL_distance(str1, str2):
    """
        Input:
            str1, str2: two strings of the same length
            representing the 
            two input structures
            in dot bracket notation e.g: '..((..))..'
        Output:
            a real or integer

    """

    tree1 = build_tree('('+str1+')')
    tree2 = build_tree('('+str2+')')

    ils1 = list_internal_leafsets(tree1)
    ils2 = list_internal_leafsets(tree2)
    
    ils1 = set([frozenset(sorted(L)) for L in ils1])
    ils2 = set([frozenset(sorted(L)) for L in ils2])

    return len(ils1.symmetric_difference(ils2))

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
            if not child[0].find('UNKNOWN') >=0:
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

    # overarching base-pair, should be only one (...(..)..(..)...)
    #                                           ^                ^
    L0 = depthk_bps(structure)
    assert(len(L0)==1)

    root = Node()
    i,j = L0[0]
    Li = Node(label=i)
    Lj = Node(label=j)
    root.add_child(Li) # other child added at the end of the function

    # rest of the leaf set
    IL = list(range(i+1,j,1))

    # second-level base-pairs (...(..)..(..)...)
    #                             ^  ^  ^  ^ 
    L1 = depthk_bps(structure, depth=1)

    # iterating over the positions within overarching bp
    pos = i+1
    while pos < j:
        # see if reached a bp of the second-level
        found = False
        for k, l in L1:
            assert(structure[k]=='(')
            assert(structure[l]==')')
            if pos==k:
                found = True
                bp_to_build = (k,l)
                break

        if found:
            k,l = bp_to_build
            # recursively building subtree
            sub_structure = list(structure)
            for index in range(len(sub_structure)):
                if index < k or index > l:
                    sub_structure[index] = '.'
            child_kl = build_tree(''.join(sub_structure))
            root.add_child(child_kl)

            # jumping
            pos = l+1
            continue

        # if none contain it, create leaf
        Lpos = Node(label=pos)
        root.add_child(Lpos)

        pos += 1

    root.add_child(Lj)
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

    # CHECK: asserting that all structures are of the same length
    lengths = set([len(val) for _, val in str_dict.items()])
    assert(len(lengths)==1)
    # END CHECK

    # building tree objects for each RNA, and getting full set of clades
    tree_dict = {}
    clade_set = set([])  

    for key, value in str_dict.items():
        value = '('+value+')' # virtual overlapping arc 
        check_well_parenthesized(value)
        tree_dict[key] = build_tree(value)
        for c in list_clades(tree_dict[key]):
            if len(c)<len(value):
                clade_set.add(frozenset(c))
        
    clade_set = list(clade_set) # order is random here
    
    trait_vector = VectorTraits(len(clade_set))

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
    '''
    Input:
        - trait_vector: 0s and 1s corresponding to the presence/absence of 
        clades of clade_set
        - clade_set: ordered list of clades
        - N: RNA length (number of positions)
    Output:
        structure
    '''

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

def best_leaflabel_solution(phylo_T, str_dict, distance=IL_distance):

    # tables
    c = {}
    distance_table = {}

    # leaf labels
    leaf_labels = []
    structure_set = set([])
    for label, structure in str_dict.items():
        if structure not in structure_set:
            structure_set.add(structure)
            leaf_labels.append(label)

    def distance_fun(s1,s2):
        try:
            return distance_table[(s1,s2)]
        except KeyError:
            distance_table[(s1,s2)] = distance(s1,s2)
            return distance_table[(s1,s2)]

    def optimal(node, leaf1, depth=0):
        """
            - node: a node (PhyloNode object)
            - leaf1: a leaf label (string)
        """
#        print(depth*'\t','optimal call on',node.label,'and',leaf1)
        # if solution cached
        if (node.label,leaf1) in c.keys():
            return c[(node.label,leaf1)]
    
        # if leaf
        if node.isleaf:
            if str_dict[node.label]==str_dict[leaf1]:
                c[(node.label, leaf1)] = 0
            else:
                c[(node.label, leaf1)] = np.inf
            return c[(node.label, leaf1)]

        c[(node.label,leaf1)] = 0
        for child in node.children:
            score = np.inf
            for leaf2 in leaf_labels:
                score = min(score, distance_fun(str_dict[leaf1], str_dict[leaf2])+optimal(child, leaf2,depth=depth+1))
            c[(node.label,leaf1)] += score

        return c[(node.label,leaf1)]

    def backtrace(node, leaf1):
        """
            - node: PhyloNode

            decides the label of the children given the ``information/constraint'' that
            node is label by leaf.
        """
        
        for child in node.children:
            # finding best leaf2 for child given leaf1 for node
            opt = np.inf
            for leaf2 in leaf_labels:
                score = distance_fun(str_dict[leaf1], str_dict[leaf2]) + c[(child.label, leaf2)]
                if score < opt:
                    best_leaf2 = leaf2
                    str_dict[child.label] = str_dict[best_leaf2]
                    opt = score

            # recursive call
            backtrace(child, best_leaf2)

    OPT = np.inf
    for leaf in leaf_labels:
        if optimal(phylo_T, leaf) < OPT:
            OPT = c[(phylo_T.label, leaf)]
            best_root_leaf = leaf
            str_dict[phylo_T.label] = str_dict[best_root_leaf]

    backtrace(phylo_T, best_root_leaf)

    return str_dict
            
    
def median_based_heuristic(phylo_T, 
                           str_dict, 
                           rounds=10, 
                           until_converge=False, 
                           median_function=C2_ILMedian,
                           distance=IL_distance,
                           best_leaf_heuristic_only=False):
    """
        - phylo_T: phylogenetic tree.
        - str_dict: dictionary label -> str. only on leaves at first, to complete
        and return
    """
    
#    random_key = np.random.choice(list(str_dict.keys()))
#    random_input_str = str_dict[random_key]
#

    str_dict = best_leaflabel_solution(phylo_T, str_dict,distance=distance)
    if best_leaf_heuristic_only:
        return str_dict

    # init: first iteration will see every node as if just changed.
    next_changed = {}
    queue = [phylo_T]
    while len(queue) > 0:
        node = queue.pop()
        next_changed[node.label] = True
        for child in node.children:
            queue.append(child)
    

    cnt = 0
    keep_going = True
    while keep_going:
        keep_going = False
        cnt += 1
        print('COUNT', cnt)

        smth_changed = False
    
        # getting info: has a node changed in last round ?
        changed = {k:v for k,v in next_changed.items()}
        print([k for k,v in changed.items() if v])
        next_changed = {k:False for k in next_changed.keys()}

        # handling the case of the root
        current_cost = 0
        for c in phylo_T.children:
            current_cost += distance(str_dict[c.label], str_dict[phylo_T.label])
        str_list = [str_dict[c.label] for c in phylo_T.children]
        median = median_function(str_list)
        new_cost = 0
        for c in phylo_T.children:
            new_cost += distance(str_dict[c.label], median)

        if new_cost < current_cost:
            str_dict[phylo_T.label] = median
            smth_changed = True

        queue = [(phylo_T, c) for c in phylo_T.children]
        while len(queue) > 0:
            parent, child = queue.pop()
            look_at_it = any([changed[parent.label]]+[changed[c.label] for c in child.children])
            print([changed[parent.label]]+[changed[c.label] for c in child.children])
            if len(child.children) > 0 and look_at_it:
                # HERE look at node and all neighbors. keep going only if one of them has ''changed''
                

                # current sum over edges around node
                current_cost = distance(str_dict[parent.label], 
                                           str_dict[child.label])
                for grand_child in child.children:
                    current_cost += distance(str_dict[child.label],
                                                str_dict[grand_child.label])

                # median computation
                neighbor_list = [parent] + child.children
                str_list = [str_dict[n.label] for n in neighbor_list]
                median = median_function(str_list)
                print('median in input structures', median in str_list)
                print('median is original structure', median==str_dict[child.label])

                # new sum over edges around node
                new_cost = distance(str_dict[parent.label], 
                                       median)
                for grand_child in child.children:
                    new_cost += distance(median,
                                                str_dict[grand_child.label])
                
                # if new cost better, change label of node to median
                print('current_cost:', current_cost)
                print('new_cost:', current_cost)
                if new_cost < current_cost:
                    str_dict[child.label] = median
                    next_changed[child.label] = True
                    # HERE, could be a step that marks this as ''changed''.
                    smth_changed = True

            
            for grand_child in child.children:
                queue.append((child, grand_child))

        # do we keep going ? count condition, and whether smth changed
        if (cnt < rounds or until_converge) and smth_changed:
            keep_going = True

    return str_dict

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


def preprocess_rfam_ncRNA_tree(ncrna_tree_line,verbose=False):
	if verbose:
		print("\n#####################################")
		print("### Correction of the ncRNA tree ####")
		print("#####################################")

	### PROBLEM WITH "(.*)" in leaf names -> cause error in ete3.Tree() reading
	# replace "(.*)" by "{.*}"
	# LEAF EXAMPLES:
	#                                    *     *
	# 66.80_S80743.1/84-123_Homo_sapiens_(human)[9606].1:0.00055
	#                                      *        *
	# 82.00_AC147908.3/69916-69833_Xenopus_(Silurana)_trop..[8364].1:0.15900
	if verbose:
		print("\n+ Replace \"(.*)\" -> \"{.*}\" :")
	i=0
	for pattern in re.compile('[^,\(\n]\(([^\)\[]{0,50}?)\)').findall(ncrna_tree_line):
		if verbose>1:
			print("\t("+pattern+") -> {"+pattern+"}")
		ncrna_tree_line=ncrna_tree_line.replace("("+pattern+")","{"+pattern+"}")
		i+=1
	if verbose:
		print("COUNT: "+str(i))


	### PROBLEM WITH "(.*[" in leaf names -> cause error in ete3.Tree() reading
	# replace "(.*[" by "{.*}["
	# LEAF EXAMPLE:
	#                                            *
	# 52.80_AJ318074.1/3315-3429_Cucumis_sativus_(cucumb..[3659].1:0.13306
	if verbose:
		print("\n+ Replace \"(.*[\" -> \"{.*}[\"")
	i=0
	for pattern in re.compile('[^,\(\n]\(([^\[]{0,50}?)\[').findall(ncrna_tree_line):
		if verbose>1:
			print("\t("+pattern+" -> {"+pattern+"}")
		ncrna_tree_line=ncrna_tree_line.replace("("+pattern,"{"+pattern+"}")
		i+=1
	if verbose:
		print("COUNT: "+str(i))

	### PROBLEM WITH ")[" in leaf names -> cause error in ete3.Tree() reading
	# replace ")[" by "}["
	# ONLY 1 LEAF EXAMPLE:
	#                                *
	# _U15164.1/177-216_spider_monkey)[36234].1:0.02519
	# replace ")[" by "}["
	if verbose:
		print("\n+ Replace \")[\" -> \"}[\"")
	ncrna_tree_line=ncrna_tree_line.replace(")[","}[")
	


	### PROBLEM WITH ":" in leaf names -> cause error in ete3.Tree() reading
	# replace ":" by "-"
	# LEAF EXAMPLE:
	#                                                  *
	# (77.20_AB035926.1/5989-6104_Escherichia_coli_O157:H7[83334].1:0.00169
	#                                                         *
	# ,231.60_FM180568.1/2999022-2999384_Escherichia_coli_O127:H..[574521].1:0.0
	if verbose:
		print("\n+ Replace \":\" -> \"-\"")
	i=0
	for pattern in re.compile('[,\(][^\[:]{0,80}?:([^\[]{0,50}?)\[').findall(ncrna_tree_line):
		if verbose>1:
			print("\t:"+pattern+" -> {"+pattern+"}")
		ncrna_tree_line=ncrna_tree_line.replace(":"+pattern,"{"+pattern+"}")
		i+=1
	if verbose:
		print("COUNT: "+str(i))


	### PROBLEM WITH ":" in leaf names whit no "[ncbi_id]"" -> cause error in ete3.Tree() reading
	# replace ":.*:" by "{.*}:"
	# LEAF EXAMPLE:
	#                                                                                   *   
	# (67.20_URS0000D6D3EB_1173780/1-51_Salmonella_enterica_subsp._diarizonae_serovar_60:r.1:0.02084
	#                                                                                *  *  
	# ,59.40_URS0000D69B8A_882884/1-50_Salmonella_enterica_subsp._arizonae_serovar_62:z4,.1:0.07437
	if verbose:
		print("\n+ Replace \":.*:\" -> \"{.*}:\"")
	i=0
	for pattern in re.compile('[,\(][^\(\):]{0,90}?:([^\(\)\[:]{0,25}?):').findall(ncrna_tree_line):
		if verbose>1:
			print(":"+pattern+" -> {"+pattern+"}")
		ncrna_tree_line=ncrna_tree_line.replace(":"+pattern,"{"+pattern+"}")
		i+=1
	if verbose:
		print("COUNT: "+str(i))

	### PROBLEM WITH ":" in leaf names -> cause error in ete3.Tree() reading
	# replace ":" by "-"
	# LEAF EXAMPLE:
	#                                                                                *   *
	# ,46.8_URS0001A238E9_1386967/1-26_Salmonella_enterica_subsp._arizonae_serovar_62:z36:-_str._RKS2983[1386967].1:0.0
	if verbose:
		print("\n+ Replace \":\" -> \"-\"")
	i=0
	for pattern in re.compile('[,\(]([^,\[:]{0,90}?:[^:,\)]{0,40}?):').findall(ncrna_tree_line):
		modif_pattern=pattern
		if ":" in pattern:
			modif_pattern=pattern.replace(':','-')
		if verbose>1:
			print(pattern+" -> "+modif_pattern)
		ncrna_tree_line=ncrna_tree_line.replace(pattern,modif_pattern)
		i+=1
	if verbose:
		print("COUNT: "+str(i))

	### PROBLEM WITH "," in leaf names -> cause error in ete3.Tree() reading
	# replace "," by "\'"
	# LEAF EXAMPLES:
	#                                                        *
	# 89.40_ACBA01000078.1/12285-12195_Enterococcus_faecium_1,..[565660].1:0.14182
	#                                                                               * *
	# 64.00_URS0000D6AD33_1173951/1-51_Salmonella_enterica_subsp._salamae_serovar_58:l,z1.1:0.02024
	if verbose:
		print("\n+ Replace \",\" -> \"\'\"")
	i=0
	for pattern in re.compile(',([^,\[:]{0,90}?,[^:]{0,25}?):').findall(ncrna_tree_line):
		modif_pattern=pattern
		if "," in pattern:
			modif_pattern=pattern.replace(',','\'')
		if verbose>1:
			print(pattern+" -> "+modif_pattern)
		ncrna_tree_line=ncrna_tree_line.replace(pattern,modif_pattern)
		i+=1
	if verbose:
		print("COUNT: "+str(i))

	### PROBLEM leaf with no RNA id and species name...
	# in all cases it seems they are multifurcated nodes
	# -> remove leaf and ","" or ")" following character 
	for pattern in re.compile('[^_]([,\(]\.\.\[[0-9]+\]\.[0-9]+\:[0-9]+\.[0-9]+[,\)])').findall(ncrna_tree_line):
		if verbose>1:
			print(pattern+" -> "+pattern[0])
		ncrna_tree_line=ncrna_tree_line.replace(pattern,pattern[0])


	# # find and remove in-parenthesis not closed in species name: "_(.*\.\."
	# for pattern in compile('(_\([^)]{0,50}?\.\.)').findall(ncrna_tree_line):
	# 	if verbose>1:
	#		print(pattern)
	# 	ncrna_tree_line=ncrna_tree_line.replace(pattern,"")

	# # A special case of in-parenthesis not closed in species name
	# # where there is no 
	# # -> Escherichia_coli_BL21(D..
	# for pattern in compile('\(([^)]{0,10}?\.\.)').findall(ncrna_tree_line):
	# 	if verbose>1:
	#		print(pattern)
	# 	ncrna_tree_line=ncrna_tree_line.replace(pattern,"{"+pattern+"}")
	# 	if verbose>1:
	#		print(ncrna_tree_line)

	# # find and remove "_(.*)" in species name
	# for pattern in compile('_\(([^)]{0,50}?)\)').findall(ncrna_tree_line):
	# 	if verbose>1:
	#		print(pattern)
	# 	ncrna_tree_line=ncrna_tree_line.replace(pattern,"{"+pattern+"}")


	### LIST OF LEAF NAMES EXAMPLES THAT HAVE TO BE MODIFIED:
	#
	#                                                        *
	# 229.00_BA000007.2/3475602-3475964_Escherichia_coli_O157:H..[386585].1:0.0

	return ncrna_tree_line

