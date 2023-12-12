from rnadist.utils import tree_tab_file_parse

def test_parse1():
    root = tree_tab_file_parse(open('test_tree.tab').readlines())

    assert(root.label=='1')
    assert([c.label for c in root.children]==['2','3'])
    assert([c.label for c in root.children[0].children]==['4','5'])
    assert([c.label for c in root.children[1].children]==[])
    assert([c.label for c in root.children[0].children[1].children]==['6'])
    assert([c.label for c in root.children[0].children[1].children[0].children]==['7','8','9'])
