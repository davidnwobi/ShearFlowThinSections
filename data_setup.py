from collections import defaultdict
import pandas as pd
from node import Node
from element import Element
from sympy import *
import warnings

import numpy as np

b, h, t, S_y, S_z, I_y, I_z, I_yz = symbols('b h t S_y S_z I_y I_z I_yz', real=True, positive=True)
dimensions = {'b': b, 'h': h, 't': t, 'S_y': S_y, 'S_z': S_z, 'I_y': I_y, 'I_z': I_z, 'I_yz': I_yz}


def extract_nodes_and_elements(nodes_path, elements_path):
    # Extract nodes and elements from excel files
    nodes_df = pd.read_excel(nodes_path, index_col=0)
    elements_df = pd.read_excel(elements_path, index_col=0)

    # Create a list of nodes and a dictionary of elements
    nodes_list = []
    for i, node_properties in nodes_df.iterrows():
        nodes_list.append(
            Node(pos=i, y=node_properties['y'] * dimensions['b'], z=node_properties['z'] * dimensions['h']))

    elements_dict = {}

    for i, element_properties in elements_df.iterrows():
        if element_properties['Ni'] > element_properties['Nj']:
            element = Element(pos=i, Node1=nodes_list[element_properties['Nj']],
                              Node2=nodes_list[element_properties['Ni']],
                              t=element_properties['t'] * dimensions['t'])
            elements_dict[(element_properties['Nj'], element_properties['Ni'])] = element
        else:
            element = Element(pos=i, Node1=nodes_list[element_properties['Ni']],
                              Node2=nodes_list[element_properties['Nj']],
                              t=element_properties['t'] * dimensions['t'])
            elements_dict[(element_properties['Ni'], element_properties['Nj'])] = element



    # Create a graph of nodes
    nodes_graph = defaultdict(list)
    for i, element in elements_df.iterrows():
        nodes_graph[element['Ni']].append((element['Nj']))
        nodes_graph[element['Nj']].append((element['Ni']))

    return nodes_df, elements_df, nodes_list, elements_dict, nodes_graph


def extract_forces(forces_path="Forces.xlsx"):
    forces_df = pd.read_excel(forces_path, index_col=0)
    force_dict = {'S_y': list(forces_df.loc['S_y']), 'S_z': list(forces_df.loc['S_z'])}
    if not force_dict['S_y'][0]:
        dimensions['S_y'] = 0
    else:
        dimensions['S_y'] = S_y
    if not force_dict['S_z'][0]:
        dimensions['S_z'] = 0
    else:
        dimensions['S_z'] = S_z
    force_dict['S_y'][1] = force_dict['S_y'][1] * dimensions['h']
    force_dict['S_z'][1] = force_dict['S_z'][1] * dimensions['b']
    return force_dict


def convert_nodes_and_elements_to_shape(nodes, elements):
    shape = pd.DataFrame(columns=['Ni', 'Nj', 'yi', 'zi', 'yj', 'zj', 't'], index=range(0))
    for i, element in elements.iterrows():
        Ni = element['Ni']
        Nj = element['Nj']
        t = element['t'] * dimensions['t']
        yi = nodes.loc[Ni]['y'] * dimensions['b']
        zi = nodes.loc[Ni]['z'] * dimensions['h']
        yj = nodes.loc[Nj]['y'] * dimensions['b']
        zj = nodes.loc[Nj]['z'] * dimensions['h']
        shape = pd.concat([shape, pd.DataFrame(
            {'Ni': [Ni], 'Nj': [Nj], 'yi': [yi], 'zi': [zi], 'yj': [yj], 'zj': [zj], 't': [t]})])
    return shape


def calculate_beam_properties(shape):
    beam_df = pd.DataFrame(shape, columns=['y', 'z'])
    beam_df['y'] = (shape['yi'] + shape['yj']) / 2
    beam_df['z'] = (shape['zi'] + shape['zj']) / 2
    beam_df['t'] = shape['t']

    beam_df_l = []
    for i in range(len(beam_df)):
        beam_df_l.append(sqrt(((shape['yj'].iloc[i] - shape['yi'].iloc[i])) ** 2 + (
            (shape['zj'].iloc[i] - shape['zi'].iloc[i])) ** 2))
    beam_df['l'] = beam_df_l
    beam_df_b = []
    for i in range(len(beam_df)):
        beam_df_b.append(abs(
            beam_df['l'].iloc[i] * (shape['yj'].iloc[i] - shape['yi'].iloc[i]) / beam_df['l'].iloc[
                i]))
    beam_df['b'] = beam_df_b

    beam_df_h = []
    for i in range(len(beam_df)):
        beam_df_h.append(abs(
            beam_df['l'].iloc[i] * (shape['zj'].iloc[i] - shape['zi'].iloc[i]) / beam_df['l'].iloc[
                i]))
    beam_df['h'] = beam_df_h

    beam_df['A'] = beam_df['t'] * beam_df['l']
    beam_df['yA'] = beam_df['y'] * beam_df['A']
    beam_df['zA'] = beam_df['z'] * beam_df['A']

    A = beam_df['A'].sum()
    dy = nsimplify(beam_df['yA'].sum() / A, tolerance=1e-10, rational=True)
    dz = nsimplify(beam_df['zA'].sum() / A, tolerance=1e-10, rational=True)

    beam_df['yi'] = beam_df['y'] - dy
    beam_df['zi'] = beam_df['z'] - dz
    beam_df['Iy'] = beam_df['t'] * (beam_df['h'] ** 2) * beam_df['l'] / 12
    beam_df['Iz'] = beam_df['t'] * (beam_df['b'] ** 2) * beam_df['l'] / 12

    beam_df_Iyz = []
    for i in range(len(beam_df)):
        tan_phi = (shape['zj'].iloc[i] - shape['zi'].iloc[i]) / (shape['yj'].iloc[i] - shape['yi'].iloc[i])
        beam_df_Iyz.append(
            beam_df['t'].iloc[i] * beam_df['b'].iloc[i] * beam_df['h'].iloc[i] * beam_df['l'].iloc[i] / 12 * sign(
                tan_phi))

    beam_df['Iyz'] = beam_df_Iyz

    beam_df['y2A'] = beam_df['yi'] ** 2 * beam_df['A']
    beam_df['z2A'] = beam_df['zi'] ** 2 * beam_df['A']
    beam_df['yzA'] = beam_df['yi'] * beam_df['zi'] * beam_df['A']

    Iy = nsimplify(beam_df['Iy'].sum() + beam_df['z2A'].sum(), tolerance=1e-10, rational=True)
    Iz = nsimplify(beam_df['Iz'].sum() + beam_df['y2A'].sum(), tolerance=1e-10, rational=True)
    Iyz = nsimplify(beam_df['Iyz'].sum() + beam_df['yzA'].sum(), tolerance=1e-10, rational=True)

    dy = simplify(dy)
    dz = simplify(dz)
    A = simplify(A)
    Iy = simplify(Iy)
    Iz = simplify(Iz)
    Iyz = simplify(Iyz)

    return dy, dz, A, Iy, Iz, Iyz


def extract_sections(sections_path="ClosedSections.xlsx"):
    sections_df = pd.read_excel(sections_path, index_col=0)
    sections = []
    try:
        Yc = (sections_df['Yc']*dimensions['b']).tolist()
        Zc = (sections_df['Zc']*dimensions['h']).tolist()
        for i, elements in sections_df.iterrows():
            sections.append(set(elements.iloc[3:].dropna().astype(int).tolist()))
        areas = (sections_df['Area'] * dimensions['b'] * dimensions['h']).tolist()
        return sections, areas, Yc, Zc
    except:
        warnings.warn("Warning: Xc and Yc not found in ClosedSections.xlsx")
        for i, elements in sections_df.iterrows():
            sections.append(set(elements.iloc[1:].dropna().astype(int).tolist()))
        areas = (sections_df['Area']*dimensions['b']*dimensions['h']).tolist()

        return sections, areas, [0]*len(sections), [0]*len(sections)



if __name__ == '__main__':
    extract_sections()
