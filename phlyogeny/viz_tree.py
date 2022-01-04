import ete3
import pandas as pd

tree = ete3.Tree('phlyogeny/augur/refined_tree.nwk', format=1)

treestyle = ete3.TreeStyle()
treestyle.branch_vertical_margin = 1

pango = pd.read_csv('phlyogeny/pangolin_assignments.tsv', sep=',').set_index('taxon')


rename_samples = pango[~pango.index.str.startswith('KB-')]['lineage'].to_dict()
focal = pango[pango.index.str.startswith('KB-')].reset_index()
focal['focal_name'] = " " + focal['scorpio_call'].str.split().str.get(0) + ":" + focal['lineage'] + " (" + focal['taxon'] + ")"
focal = focal.set_index('taxon')
focal = focal['focal_name'].to_dict()

for node in tree:
    if node.name == "MN908947.3":
        node.name = " Ancestral SARS-CoV-2"
    elif node.name in rename_samples:
        node.name = rename_samples[node.name]
    elif node.name in focal:
        node.name = focal[node.name]

nodes_to_keep = ["B.1.1.7", "Q.8", "Q.4", "Q.3", "B.1.351", "P.1", "P.1.14", "P.1.17", "BA.1", "BA.2", "B.1.617.2", "AY.25", "AY.74", "AY.27", "AY.103", "AY.93"]
nodes_to_keep.append("Ancestral SARS-CoV-2")
#for value in focal.values():
#    nodes_to_keep.append(value)

ancestral_color = "#4a5fa7"
beta_color = "#bf6aa7"
delta_color = "#f49737"
omicron_color = "#953a92"


alpha_nodes = []
beta_nodes = []
gamma_nodes = []
delta_nodes = []
omicron_nodes = []
ancestral_nodes = []

# color internal leaves grey
for node in tree.traverse():
    if not node.is_leaf():
        nstyle = ete3.NodeStyle()
        nstyle['fgcolor'] = 'BCBCBC'
        nstyle['size'] = 2
        node.set_style(nstyle)


for node in tree.iter_leaves():
    nstyle = ete3.NodeStyle()
    if node.name in ["B.1.1.7", "Q.8", "Q.4", "Q.3"]:
        alpha_nodes.append(node)
        nstyle['fgcolor'] = 'green'
        nstyle['size'] = 5
        if node.name == 'B.1.1.7':
            node.name = "Alpha:B.1.1.7"
    elif node.name.startswith('AY.') or node.name == 'B.1.617.2' or node.name == " Delta:B.1.617.2 (KB-D-002)":
        delta_nodes.append(node)
        nstyle['fgcolor'] = delta_color
        if node.name in nodes_to_keep:
            nstyle['size'] = 5
        elif node.name == " Delta:B.1.617.2 (KB-D-002)":
            nstyle['size'] = 15
        else:
            node.name = ''
            nstyle['size'] = 3
            nstyle["hz_line_color"] = "#cccccc"
            nstyle["vt_line_color"] = "#cccccc"
    elif node.name in ["P.1", "P.1.14", "P.1.17"]:
        gamma_nodes.append(node)
        nstyle['fgcolor'] = 'brown'
        nstyle['size'] = 5
        if node.name == "P.1.14":
            node.name = "Gamma:P.1.14"
    elif node.name in ["BA.1", "BA.2", ' Omicron:BA.1 (KB-D-003)']:
        omicron_nodes.append(node)
        nstyle['fgcolor'] = omicron_color
        nstyle['size'] = 5
        if node.name == ' Omicron:BA.1 (KB-D-003)':
            nstyle['size'] = 15
    elif node.name in ["B.1.351", " Beta:B.1.351 (KB-B-001)"]:
        beta_nodes.append(node)
        nstyle['fgcolor'] = beta_color
        nstyle['size'] = 5
        if node.name == " Beta:B.1.351 (KB-B-001)":
            nstyle['size'] = 15
    elif node.name in [" Ancestral SARS-CoV-2"]:
        ancestral_nodes.append(node)
        nstyle['fgcolor'] = ancestral_color
        nstyle['size'] = 5
        if node.name == ' Ancestral SARS-CoV-2':
            nstyle['size'] = 15
    else:
        node.name = ''
        nstyle['size'] = 2
        nstyle["hz_line_color"] = "#cccccc"
        nstyle["vt_line_color"] = "#cccccc"
        nstyle["fgcolor"] = "#cccccc"
    node.set_style(nstyle)

#
#for node in tree.get_common_ancestor(alpha_nodes).traverse():
#    if not node.is_leaf():
#        nstyle = ete3.NodeStyle()
#        nstyle['fgcolor'] = 'green'
#        node.set_style(nstyle)
#
#
#for node in tree.get_common_ancestor(beta_nodes).traverse():
#    if not node.is_leaf():
#        nstyle = ete3.NodeStyle()
#        nstyle['fgcolor'] = beta_color
#        node.set_style(nstyle)
#
#for node in tree.get_common_ancestor(gamma_nodes).traverse():
#    if not node.is_leaf():
#        nstyle = ete3.NodeStyle()
#        nstyle['fgcolor'] = 'brown'
#        node.set_style(nstyle)
#
#for node in tree.get_common_ancestor(delta_nodes).traverse():
#    if not node.is_leaf():
#        nstyle = ete3.NodeStyle()
#        nstyle['fgcolor'] = delta_color
#        node.set_style(nstyle)
#
#for node in tree.get_common_ancestor(omicron_nodes).traverse():
#    if not node.is_leaf():
#        nstyle = ete3.NodeStyle()
#        nstyle['fgcolor'] = omicron_color
#        node.set_style(nstyle)
#
#for node in tree.get_common_ancestor(ancestral_nodes).traverse():
#    if not node.is_leaf():
#        nstyle = ete3.NodeStyle()
#        nstyle['fgcolor'] = ancestral_color
#        node.set_style(nstyle)
#


##pango.se
tree.render('tree.svg', tree_style=treestyle, dpi=300)
