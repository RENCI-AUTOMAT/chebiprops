import wget
import requests
import jsonlines
from gzip import GzipFile
from collections import defaultdict

def pull_file(fname):
    fdir='https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/'
    fname = f'{fdir}/{fname}'
    wget.download(fname)

def read_names():
    names = {}
    with GzipFile('compounds.tsv.gz','rb') as inf:
        for bytesline in inf:
            line = bytesline.decode('utf-8')
            x = line.strip().split('\t')
            cid = x[2]
            t = x[5]
            names[cid] = t
    return names

def update_ancestors(ancestors,parent,iia):
    kids = iia[parent]
    for kid in kids:
        ancestors[kid].append(parent)
        ancestors[kid] += ancestors[parent]
        update_ancestors(ancestors,kid,iia)

def get_ancestors(iia):
    ancestors = defaultdict(list)
    role = 'CHEBI:50906'
    update_ancestors(ancestors,role,iia)
    print(ancestors['CHEBI:50904'])
    return ancestors

def read_roles():
    """This format is not completely obvious, but the triple is (FINAL_ID)-[type]->(INIT_ID)."""
    roles = defaultdict(list)
    invert_is_a = defaultdict(list)
    with open('relation.tsv','r') as inf:
        for line in inf:
            x = line.strip().split('\t')
            if x[1] == 'has_role':
                roles[f'CHEBI:{x[3]}'].append(f'CHEBI:{x[2]}')
            elif x[1] == 'is_a':
                child = f'CHEBI:{x[3]}'
                parent = f'CHEBI:{x[2]}'
                invert_is_a[parent].append(child)
    #Now include parents
    ancestors = get_ancestors(invert_is_a)
    for node,noderoles in roles.items():
        if node == 'CHEBI:64663':
            print('hi')
        restroles= []
        for role in noderoles:
            moreroles=ancestors[role]
            restroles += moreroles
        roles[node] += restroles
    return roles

def chunk(l,n):
    for i in range(0,len(l),n):
        yield l[i:i + n]

def normalize_all(statuses):
    nodes = {}
    for idlist in chunk(list(statuses.keys()), 100):
        nodes.update(normalize(idlist))
    return nodes

def normalize(il):
    nodenorm_url = 'https://nodenormalization-sri.renci.org/get_normalized_nodes'
    input = { 'curies': il }
    results = requests.post(nodenorm_url,json = input)
    return results.json()

def fixname(n):
    return f'CHEBI_ROLE:{"_".join(n.split())}'

def transform(nodes,roles,names):
    transformed = []
    for node_id, node in nodes.items():
        if node is None:
            continue
        outnode = {'id': node['id']['identifier'], 'name': node['id'].get('label',''),
                   'equivalent_identifiers': [eid['identifier'] for eid in node['equivalent_identifiers']],
                   'category': node['type']}
        rolenames = [ fixname(names[x]) for x in roles[node_id] ]
        for rn in rolenames:
            if rn in ['CHEBI_ROLE:role', 'CHEBI_ROLE:biological_role', 'CHEBI_ROLE:chemical_role', 'CHEBI_ROLE:application']:
                continue
            outnode[rn] = True
        transformed.append(outnode)
    return transformed

def write(nodes):
    with jsonlines.open('ChemProperties.jsonl','w') as outf:
        for node in nodes:
            outf.write(node)

def go():
    #pull_file('compounds.tsv.gz')
    #pull_file('relation.tsv')
    names = read_names()
    roles = read_roles()
    print(len(roles))
    n = normalize_all(roles)
    print(len(n))
    t = transform(n,roles,names)
    print(len(t))
    write(t)

if __name__ == '__main__':
    go()