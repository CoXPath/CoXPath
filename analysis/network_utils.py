# -*- coding: utf-8 -*-

"""Helper functions when working with networks in CoXPath"""

import os
from collections import Counter
from itertools import chain, product
from typing import Iterable, List, Dict, Union, Tuple, Set

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import pyobo
from matplotlib_venn import venn2

homedir = os.path.expanduser("~")
network_cache_dir = os.path.join(homedir, ".coexp_nets")
os.makedirs(network_cache_dir, exist_ok=True)

doid_mappings = pyobo.get_id_name_mapping('doid')
hgnc_mappings = pyobo.get_id_name_mapping('hgnc')
hgnc_id_to_ncbigene_id = pyobo.get_filtered_xrefs('hgnc', 'ncbigene')
entrez_mappings = {
    ncbigene_id: hgnc_mappings[hgnc_id]
    for hgnc_id, ncbigene_id in hgnc_id_to_ncbigene_id.items()
}

class InputError(Exception):
    """Exception raised for errors in the input.

    Attributes
    ----------
    message : str
        explanation of the error
    """

    def __init__(self, message):
        self.message = message


def edge_file_path(data_dir: str, doid: str) -> str:
    """Return the path for the edge list file for a given doid

    Parameters
    ----------
    data_dir : str
        path to the directory which contains co-expression network edges.
    doid : str
        disease DOID

    Returns
    -------
        path to edge file
    """
    return os.path.join(data_dir, doid, "coexp_network_edges.tsv")


def create_network_from_edge_file(file: str, doid: str) -> nx.Graph:
    """Extract the necessary columns to create a weighted undirected graph

    Parameters
    ----------
    file : str
        path to edge list life
    doid : str
        disease DOID

    Returns
    -------
        network as nx.Graph object
    """
    try:  # see if it's already cached
        G = nx.read_gpickle(os.path.join(network_cache_dir, f"{doid}.gpickle"))
    except FileNotFoundError:  # if not, create and cache the graph for next time
        edge_df = pd.read_csv(file, sep="\t", usecols=["from", "to", "weight"])
        G = nx.from_pandas_edgelist(edge_df, source="from", target="to", edge_attr="weight")
        nx.write_gpickle(G, os.path.join(network_cache_dir, f"{doid}.gpickle"))
    return G


def counter(nets: Iterable, comparison: str) -> Counter:
    """Create a counter objects of nodes/edges/genes from a collection of Networks.

    Should not be called by user directly- only by other methods.
    """
    if comparison == "nodes":
        comp = (net.nodes() for net in nets)
    elif comparison == "edges":
        comp = (net.edges() for net in nets)
    elif comparison == "genes":
        comp = nets
    else:
        raise InputError("Can only compare 'nodes' , 'edges', or 'genes'")
    counter_obj = Counter(chain.from_iterable(comp))
    return counter_obj


def most_common(nets: Iterable, comparison: str) -> List[Tuple[Union[str, tuple], Tuple[str, str]]]:
    """Find most common nodes or edges in a collection of networks

    Paramaters
    ----------
    nets : Iterable
        a collection containing the Networkx objects to analyze or of lists of gene names
    comparison : str
        type of object being compared. only accepts 'nodes', 'edges', or 'genes'

    Returns
    -------
    List[Tuple[Union[str,tuple],Tuple[str,str]]]
        list of (comparison_object, (frequency_count, percent)) tuples in order of most common
        comparison_object being nodes(str), edges(tuple of str), or genes(str)
        note - will contain all nodes/edges/genes in this list, check the counts for the info you need
    """
    counter_obj = counter(nets, comparison)
    return [(obj, (f"{count}/{len(nets)}", f"{count / len(nets):.2%}")) for obj, count in counter_obj.most_common()]


def least_common(nets: Iterable, comparison: str) -> List[Tuple[Union[str, tuple], Tuple[str, str]]]:
    """Find least common nodes or edges in a collection of networks

    Paramaters
    ----------
    nets : Iterable
        a collection containing the Networkx objects to analyze or of lists of gene names
    comparison : str
        comparison type. only accepts 'nodes', 'edges', or 'genes'

    Returns
    -------
    List[Tuple[Union[str,tuple],Tuple[str,str]]]
        list of (comparison_object, (frequency_count, percent)) tuples in order of least common
        comparison_object being nodes(str), edges(tuple of str), or genes(str)
        note - will contain all nodes/edges/genes in this list, check the counts for the info you need

    """
    counter_obj = counter(nets, comparison)
    return [(obj, (f"{count}/{len(nets)}", f"{count / len(nets):.2%}")) for obj, count in
            counter_obj.most_common()[::-1]]


def count_in_freq(freq_tuple: Tuple[str, str]) -> int:
    """Return the count from the statistic half of the tuple returned from most_common or least_common

    Parameters
    ----------
    freq_tuple : Tuple[str, str]
        tuple containing stats given as (count/total, percent) from most_common or least_common

    Returns
    -------
        the count from the statistic as an int

    """
    return int(freq_tuple[0].split('/')[0])


def above_cutoff(gene_freq_tup_list: List[Tuple[Union[str, tuple], Tuple[str, str]]], cutoff: int) -> List[str]:
    """Return the genes/edges that are are in at least the given cutoff's networks

    Parameters
    ----------
    gene_freq_tup_list : List[Tuple[Union[str, tuple], Tuple[str, str]]]
        list of (comparison_object, (frequency_count, percent)) tuples in order of most common
        should be return from most_common()
    cutoff : int
        number to be used as minimum for how many networks the object must be present in to be returned

    Returns
    -------
        list of objects that were in at least as many networks as the cutoff given
    """
    above = []
    for gene, freq in gene_freq_tup_list:
        if count_in_freq(freq) >= cutoff:
            above.append(gene)
        else:
            break  # since it's ordered, no need wasting time checking the rest
    return above


def below_cutoff(gene_freq_tup_list: List[Tuple[Union[str, tuple], Tuple[str, str]]], cutoff: int) -> List[str]:
    """Return the genes/edges that are are in less than the given cutoff's networks

    Parameters
    ----------
    gene_freq_tup_list : List[Tuple[Union[str, tuple], Tuple[str, str]]]
        list of (comparison_object, (frequency_count, percent)) tuples in order of least common
        should be return from least_common()
    cutoff : int
        number to be used as maximum for how many networks the object must be present in to be returned

    Returns
    -------
        list of objects that were in at most as many networks as the cutoff given
    """
    below = []
    for gene, freq in gene_freq_tup_list:
        if count_in_freq(freq) <= cutoff:
            below.append(gene)
        else:
            break  # since it's ordered, no need wasting time checking the rest
    return below


def identify_hub_genes(net: nx.Graph) -> List[Tuple[str, int]]:
    """Identify the nodes with the highest connections

    Parameters
    ----------
    net : nx.Graph
        network object

    Returns
    -------
        list of all nodes and their degrees, ordered by highest number of connections first

    """
    degree_sequence = sorted(net.degree, key=lambda x: x[1], reverse=True)
    return degree_sequence


def load_interactome(path_to_interactome_file: str) -> nx.DiGraph:
    """Load the interactome located at the given file location

    Parameters
    ----------
    path_to_interactome_file: str
        path to the location of the interactome file

    Returns
    -------
        interactome given as nextworkx directed graph object

    """
    interactome_df = pd.read_csv(path_to_interactome_file, sep="\t")
    G = nx.DiGraph()
    for index, row in interactome_df.iterrows():
        G.add_edge(hgnc_mappings[str(row.source)], hgnc_mappings[str(row.target)])
    return G


def get_genes_in_overlap(genes, overlap_list) -> Set[str]:
    """Select members from a gene list which are also in a given overlap(second list of genes)

        Parameters
        ----------
        genes : list
            ordered list of genes from which to select from
        overlap_list : list
            list of genes from which selected genes must be in
        cutoff

        Returns
        -------
            list of genes that are also present in the overlap list
        """
    top_genes = set()
    for gene in genes:
        if gene in overlap_list:
            top_genes.add(gene)
    return top_genes


def get_top_x_genes_in_overlap(genes: list, overlap_list: list, topx: int) -> Set[str]:
    """Select a certain number of genes from an ordered list which are also in a given overlap(second list of genes)

    Parameters
    ----------
    genes : list
        list of (genes, stats) from which to take the top x genes from
    overlap_list : list
        list of genes
    topx : int
        number of desired genes in resulting list

    Returns
    -------
        list of genes of length (topx) that are also present in the overlap list
    """
    index = 0
    top_genes = set()
    while len(top_genes) < topx:
        gene, _ = genes[index]
        if gene in overlap_list:
            top_genes.add(gene)
        index += 1
    return top_genes


def venn_d(left: set, right: set, left_name: str, right_name: str, title: str = None,
           output_filename: str = None) -> venn2:
    """Analyze node overlap in two NetworkX objects

    Parameters
    ----------
    left : set
        set of objects for left side of venn diagram
    right : set
        set of objects for right side of venn diagram
    left_name : str
        name of collection on left side of venn diagram
    right_name : str
        name of collection on right side of venn diagram
    title : str
        title of venn diagram, optional
    output_filename : str
        path to desired location of saved figure, optional

    Returns
    -------
        venn diagram object

    """
    v = venn2([left, right], (left_name, right_name))
    v.get_patch_by_id('10').set_color('#5cb85c')
    v.get_patch_by_id('01').set_color('#df3f18')
    v.get_patch_by_id('10').set_edgecolor('none')
    v.get_patch_by_id('01').set_edgecolor('none')
    v.get_patch_by_id('10').set_alpha(.3)
    v.get_patch_by_id('01').set_alpha(.3)
    if v.get_patch_by_id('11'):
        v.get_patch_by_id('11').set_color('#f3ac1f')
        v.get_patch_by_id('11').set_edgecolor('none')
        v.get_patch_by_id('11').set_alpha(.3)
    for text in v.set_labels:
        text.set_fontsize("small")
    for x in range(len(v.subset_labels)):
        if v.subset_labels[x] is not None:
            v.subset_labels[x].set_fontsize("small")
    if title:
        plt.title(title)
    if output_filename:
        plt.savefig(output_filename)
    return v


def node_intersection_from_doid_list(doid_list: List[str], network_d: Dict[str, nx.Graph]) -> Set[str]:
    """Find node intersection in the networks of the given doid list

    Parameters
    ----------
    doid_list: List[str]
        list of doids
    network_d: Dict[str, nx.Graph]
        dictionary mapping doids to their network

    Returns
    -------
        set of nodes
    """
    nodes = [set(network_d[doid].nodes()) for doid in doid_list]
    return nodes[0].intersection(*nodes[1:])


def remove_stats(obj_stat_list: List[Tuple[Union[str, tuple], Tuple[str, str]]]) -> List[Tuple[Union[str, tuple]]]:
    """Given the return from most_common or least_common, remove the stats and just return a list of the nodes/edges

    Parameters
    ----------
    obj_stat_list:
        list of (comparison_object, (frequency_count, percent)) tuples in which you only want the comparison objects

    Returns
    -------
        list of objects (can either be edges or nodes)
    """
    return [obj for obj, stat in obj_stat_list]


def percent_edges_in_other_edges(edge_list: List[Tuple[str, str]], edge_list2: List[Tuple[str, str]]) -> float:
    """Calculate what proportion of the first list of edges can be found in the second list.

    checks edges in both node orders to account for directed edges

    Parameters
    ----------
    edge_list:  List[Tuple[str, str]]
        list of query edges
    edge_list2:  List[Tuple[str, str]]
        list of edges to search in

    Returns
    -------
        proportion of first list that is in second list, accounting for order differences
    """
    count = 0
    for edge in edge_list:
        if (edge in edge_list2) or (edge[::-1] in edge_list2):
            count += 1
    return count / len(edge_list)


def top_edges_subgraph(num_edges: int, sorted_edges: List[Tuple[str, str]]) -> nx.Graph:
    """Given an ordered list(by weight) of edges from a network, create a subgraph containing only a given number of edges

    Parameters
    ----------
    num_edges : int
        desired number of edges to be in resulting subgraph
    sorted_edges: List[Tuple[str, str]]
        list of sorted edges by weight

    Returns
    -------
        subgraph networkx object

    """
    subG = nx.Graph()
    subG.add_edges_from(sorted_edges[:num_edges])
    return subG


def pathway_similarity(gene_list_in_pathway: List[str], network: nx.Graph()) -> float:
    """Calculate pathway similarity.

    Given a list of genes in a pathway, looks for edges between these genes in the given network

    Parameters
    ----------
    gene_list_in_pathway: List[str]
        list of genes in a pathway
    network: nx.Graph
        network to search for pathway similarity

    Returns
    -------
        similarity value
    """
    edges_found = 0
    i = 0
    for gene1, gene2 in product(gene_list_in_pathway, repeat=2):
        i += 1
        if network.has_edge(gene1, gene2):
            edges_found += 1
    return edges_found / i


def pathway_similarity_via_interactome(gene_list_in_pathway: List[str], network: nx.Graph(),
                                       interactome: nx.DiGraph()) -> float:
    """Calculate pathway similarity.

    Given a list of genes in a pathway, looks for edges between these genes that are also in the interactome in the given network

    Parameters
    ----------
    gene_list_in_pathway: List[str]
        list of genes in a pathway
    network: nx.Graph
        network to search for pathway similarity
    interactome: nx.DiGraph

    Returns
    -------
        similarity value
    """
    interactome_subgraph_with_pathway = interactome.subgraph(
        [node for node in gene_list_in_pathway if node in interactome.nodes()])
    edges_found = 0
    if not interactome_subgraph_with_pathway.edges():
        return 0
    for n1, n2 in interactome_subgraph_with_pathway.edges():
        if network.has_edge(n1, n2):
            edges_found += 1
    return edges_found / len(interactome_subgraph_with_pathway.edges())


def largest_subgraph(G1: nx.Graph, G2: nx.Graph) -> nx.Graph:
    """Return largest connected subgraph that is overlapping fully between two graphs

    Parameters
    ----------
    G1 : nx.Graph
        first network
    G2 : nx.Graph
        second network

    Returns
    -------
        largest common subgraph
    """
    intersectionG = nx.intersection(G1.subgraph(G2.nodes), G2.subgraph(G1.nodes))  # still contains many "empty" nodes
    subG = nx.Graph()
    subG.add_edges_from(intersectionG.edges)

    # not necessarily connected yet
    if subG:
        max_comp = max(nx.connected_components(subG), key=len)
        subG = subG.subgraph(max_comp)
    return subG


def largest_subgraph_from_list(nets: List[nx.Graph]) -> nx.Graph:
    """Return largest connected subgraph that is overlapping fully between many graphs

    Parameters
    ----------
    nets : List[nx.Graph]
        list of networks

    Returns
    -------
        largest common subgraph
    """
    for i, net in enumerate(nets):
        if i == 0:
            subg = net
            continue
        subg = largest_subgraph(subg, net)
    return subg


def which_nets_have(item: Union[str, tuple], item_type: str, network_dict: Dict[str, nx.Graph]) -> List[str]:
    """Identify which network contains a given node/edge

    Parameters
    ----------
    item: Union[str,tuple]
        node(str) or edge(tuple) that you want to check which networks have
    item_type : str
        whether to search for a node or edge. only accepts 'node' or 'edge'
    network_dict

    Returns
    -------
        list of doids of the networks that have the query

    """
    if item_type == "node":
        search = {doid: net.nodes() for doid, net in network_dict.items()}
    elif item_type == "edge":
        search = {doid: net.edges() for doid, net in network_dict.items()}
    else:
        raise InputError("Can only find 'node' and 'edge'")
    doids_that_have_item = []
    for doid, group in search.items():
        if item in group:
            doids_that_have_item.append(doid)
    return doids_that_have_item


def ensp_to_hgnc_mappings(URL = "https://stringdb-static.org/download/protein.info.v11.0/9606.protein.info.v11.0.txt.gz"):
    """Map Ensembl protein ids to HGNC symbols
    
    Parameters
    ----------
    URL: str
        location of the chosen version's STRING protein info file. taken from: https://string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Homo+sapiens
        version 11 human protein info is given as the default, optional

    Returns
    -------
        dictionary mapping ensembl ids to hgnc gene symbols
    
    """
    ensembl_hgnc_mappings_df = pd.read_table(URL, compression ="gzip")

    ensp_gene_mappings = {
        row["protein_external_id"] : row["preferred_name"]
        for i, row in ensembl_hgnc_mappings_df.iterrows()
    }
    return ensp_gene_mappings


def load_STRING(URL = "https://stringdb-static.org/download/protein.physical.links.v11.0/9606.protein.physical.links.v11.0.txt.gz"):
    """Load the STRING PPI network located at the given URL

    Parameters
    ----------
    URL: str
        location of the chosen version's STRING PPI physical links file. taken from: https://string-db.org/cgi/download?sessionId=%24input-%3E%7BsessionId%7D&species_text=Homo+sapiens
        version 11 human links are given as the default, optional

    Returns
    -------
        PPI network given as nextworkx directed graph object

    """
    STRING_edges = pd.read_table(URL, compression ="gzip", sep=" ")
    STRING = nx.Graph()
    gene_mappings = ensp_to_hgnc_mappings()
    for i, row in STRING_edges.iterrows():
        source = gene_mappings[row["protein1"]]
        target = gene_mappings[row["protein2"]]
        STRING.add_edge(source, target, weight = row["combined_score"])
    return STRING


def load_HIPPIE(URL = "http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt"):
    """Load the HIPPIE PPI network located at the given URL

    Parameters
    ----------
    URL: str
        location of the chosen version's HIPPIE's PPI tab format file. taken from: http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php
        the link automatically uses the current version (2.2  at the time this was conducted) as the default, optional

    Returns
    -------
        PPI network given as nextworkx directed graph object

    """
    HIPPIE_interactions = pd.read_table(URL, header=None)
    HIPPIE= nx.Graph()
    skipped = 0
    for i, row in HIPPIE_interactions.iterrows():
        try:
            source = entrez_mappings[str(row[1])]
            target = entrez_mappings[str(row[3])]
        except KeyError:
            skipped += 1
            continue
        HIPPIE.add_edge(source, target, weight = row[4])
    print(f"{skipped} edges skipped due to discontinued gene ids")
    return HIPPIE
