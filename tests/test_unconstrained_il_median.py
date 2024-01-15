from rnadist.indset_intgraph import get_max_interval_indset

def test_empty_input():

    output = get_max_interval_indset([])

def test_empty_output():

    output = get_max_interval_indset([(35, 36, -3)])
