import unittest
from go import GOMap

obo_path = r".\go.obo"


class TestGoMapMethods(unittest.TestCase):
    def test_parser(self):
        g_map = GOMap()
        g_map.parse_obo(obo_path)

    def test_getting_go_term(self):
        g_map = GOMap()
        g_map.parse_obo(obo_path)
        id = "GO:0000011"
        g = g_map.get_go_term(id, "biological_process")
        self.assertEqual(g.id, id)
        g = g_map.get_go_term("GO:0000084", "biological_process")
        self.assertEqual(2, len(g.parents))

    def test_getting_all_parents(self):
        g_map = GOMap()
        g_map.parse_obo(obo_path)
        id = ["GO:0000011", "GO:0000022"]
        results = g_map.get_all_parents(id, "biological_process")
        print(results)

    def test_get_shortest_paths(self):
        g_map = GOMap()
        g_map.parse_obo(obo_path)
        id = ["GO:0000011", "GO:0000022"]
        common, data, lowest_depth = g_map.get_shortest_path(id, "biological_process")
        data.to_csv("test.csv", index=False)
        print(common, data)
