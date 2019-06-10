import re
import pandas as pd
from copy import deepcopy

is_a_pattern = re.compile("(GO:\d+)")
quote_pattern = re.compile("\"(.+)\"")
base_path = {"GO:0008150", "GO:0005575", "GO:0003674"}


class GOTerm:
    def __init__(self, id=None, name=None, ns=None, definition=None):
        self.id = id
        self.name = name
        self.ns = ns
        self.definition = definition
        self.parents = []
        self.children = []

    def get_all_parents(self, go_map, pathway, current_path=None, current_depth=0):
        if not current_path:
            current_path = {}
        current_path[current_depth] = ""
        for p in self.parents:
            parent = go_map.get_go_term(p, pathway)
            current_path[current_depth] = parent.id
            if p in base_path:
                yield current_path, current_depth
            else:
                yield from parent.get_all_parents(go_map, pathway, deepcopy(current_path), current_depth=current_depth+1)


def get_go_info(obo, object_dict):
    go_term = GOTerm()
    for line in obo:
        line = line.strip()
        if line.startswith("id:"):
            go_term.id = line[4:]
        elif line.startswith("name:"):
            go_term.name = line[6:]
        elif line.startswith("namespace:"):
            if go_term.id in object_dict[line[11:]]:
                go_term = object_dict[line[11:]][go_term.id]
            go_term.ns = line[11:]
        elif line.startswith("alt_id:"):
            object_dict[go_term.ns][line[8:]] = go_term
        elif line.startswith("def:"):
            s = quote_pattern.search(line)
            if s:
                go_term.definition = s.group(1)
        elif line.startswith("is_a:"):
            s = is_a_pattern.search(line)
            if s:
                go_term.parents.append(s.group(1))
                if s.group(1) not in object_dict[go_term.ns]:
                    object_dict[go_term.ns][s.group(1)] = GOTerm(id=s.group(1))
                object_dict[go_term.ns][s.group(1)].children.append(go_term.id)
        elif line == "":
            if go_term.id not in object_dict[go_term.ns]:
                object_dict[go_term.ns][go_term.id] = go_term
            break


class GOMap:
    def __init__(self):
        self.object_dict = {"biological_process": {}, "molecular_function": {}, "cellular_component": {}}

    def parse_obo(self, obo_path):
        with open(obo_path, "rt") as obo:
            for line in obo:
                line = line.strip()
                if line.startswith("[Term]"):
                    get_go_info(obo, self.object_dict)

    def get_go_term(self, go_id: str, pathway: str) -> GOTerm:
        return self.object_dict[pathway][go_id]

    def get_all_parents(self, go_ids, pathway):
        results = []
        lowest_depth = 0
        for i in go_ids:
            path_count = 0
            for path, depth in self.get_go_term(i, pathway).get_all_parents(self, pathway):
                path_count += 1
                if depth > lowest_depth:
                    lowest_depth = depth
                row = {"query": i, "path_number": path_count}
                for r in range(0, depth+1, 1):
                    row[r] = path[depth - r]
                results.append(row)

        return pd.DataFrame(results), lowest_depth

    def get_shortest_path(self, go_ids, pathway):
        paths, lowest_depth = self.get_all_parents(go_ids, pathway)
        previous_term = ""
        previous_depth = 0
        for d in range(0, lowest_depth+1, 1):
            unique_term = {}
            for i, g in paths.groupby("query"):
                unique_term[i] = g[d].unique()

            for t in unique_term[i]:
                occur = 0
                for c in unique_term:
                    if c != t:
                        if t in unique_term[c]:
                            occur += 1
                if occur == len(go_ids):
                    previous_term = t
                    previous_depth = d
            if d != previous_depth:
                return previous_term, paths, previous_depth







