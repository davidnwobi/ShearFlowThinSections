import pandas as pd
import numpy as np


def extract_nodes_and_elements(nodes_path="Nodes.xlsx", elements_path="Elements.xlsx"):
    nodes_df = pd.read_excel(nodes_path, index_col=0)
    elements_df = pd.read_excel(elements_path, index_col=0)

    return nodes_df, elements_df

def convert_nodes_and_elements_to_shape(nodes, elements):
    shape = pd.DataFrame(columns=['Ni', 'Nj','yi', 'zi', 'yj', 'zj', 't'], index=range(0))
    for i, element in elements.iterrows():
        Ni = element['Ni']
        Nj = element['Nj']
        t = element['t']
        yi = nodes.loc[Ni]['y']
        zi = nodes.loc[Ni]['z']
        yj = nodes.loc[Nj]['y']
        zj = nodes.loc[Nj]['z']
        shape = pd.concat([shape, pd.DataFrame({'Ni': [Ni], 'Nj': [Nj], 'yi': [yi], 'zi': [zi], 'yj': [yj], 'zj': [zj], 't': [t]})])
    return shape


def calculate_beam_properties(file_path="shape.xlsx"):
    df = pd.read_excel("shape.xlsx", index_col=0)
    beam_df = pd.DataFrame(df, columns=['y', 'z'])
    beam_df['y'] = (df['yi'] + df['yj']) / 2
    beam_df['z'] = (df['zi'] + df['zj']) / 2
    beam_df['t'] = df['t']
    beam_df['phi'] = np.arctan2(df['zj'] - df['zi'], df['yj'] - df['yi'])
    beam_df['l'] = np.sqrt((df['yj'] - df['yi']) ** 2 + (df['zj'] - df['zi']) ** 2)
    beam_df['b'] = np.abs(beam_df['l'] * np.cos(beam_df['phi']))
    beam_df['h'] = np.abs(beam_df['l'] * np.sin(beam_df['phi']))
    beam_df['A'] = beam_df['t'] * beam_df['l']
    beam_df['yA'] = beam_df['y'] * beam_df['A']
    beam_df['zA'] = beam_df['z'] * beam_df['A']

    A = beam_df['A'].sum()
    dy = beam_df['yA'].sum() / A
    dz = beam_df['zA'].sum() / A

    beam_df['yi'] = beam_df['y'] - dy
    beam_df['zi'] = beam_df['z'] - dz
    beam_df['Iy'] = beam_df['t'] * (beam_df['h'] ** 2) * beam_df['l'] / 12
    beam_df['Iz'] = beam_df['t'] * (beam_df['b'] ** 2) * beam_df['l'] / 12
    beam_df['Iyz'] = beam_df['t'] * beam_df['b'] * beam_df['h'] * beam_df['l'] / 12 * np.sign(np.tan(beam_df['phi']))
    beam_df['y2A'] = beam_df['yi'] ** 2 * beam_df['A']
    beam_df['z2A'] = beam_df['zi'] ** 2 * beam_df['A']
    beam_df['yzA'] = beam_df['yi'] * beam_df['zi'] * beam_df['A']
    Iy = beam_df['Iy'].sum() + beam_df['z2A'].sum()
    Iz = beam_df['Iz'].sum() + beam_df['y2A'].sum()
    Iyz = beam_df['Iyz'].sum() + beam_df['yzA'].sum()

    return dy, dz, A, Iy, Iz, Iyz
