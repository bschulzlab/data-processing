from go import GOMap

# ids os the goids to find the shortest common term
# pathway is the name of the pathway biological_process, molecular_function, or cellular_component
ids = ["GO:0000011", "GO:0000022"]
pathway = "biological_process"


if __name__ == "__main__":
    obo_path = r".\go.obo"
    g_map = GOMap()
    g_map.parse_obo(obo_path)

    common, data, lowest_depth = g_map.get_shortest_path(ids, pathway)
    data.to_csv("test.csv", index=False)
    print("Closest common GOTerm:{} at depth{}".format(common, lowest_depth))




