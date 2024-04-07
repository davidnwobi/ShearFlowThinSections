import pandas as pd
from ..utils.data_setup import extract_nodes_and_elements, convert_nodes_and_elements_to_shape, calculate_beam_properties, \
    extract_forces, extract_sections, dimensions


class Shape:
    def __init__(self, nodes="Nodes.xlsx", elements="Elements.xlsx", forces="Forces.xlsx", sections="ClosedSections"
                                                                                                    ".xlsx",
                 shear_center_ref="Shear Center Ref.xlsx"):
        nodes_df, elements_df, self.nodes, self.elements, self.node_graph = extract_nodes_and_elements(nodes,
                                                                                                       elements)
        shape = convert_nodes_and_elements_to_shape(nodes_df, elements_df)
        self.force_dict = extract_forces(forces)
        self.dy, self.dz, self.Ac, self.Iy, self.Iz, self.Iyz = calculate_beam_properties(shape)
        self.ey, self.ez = 0, 0
        self.sections, self.areas, self.section_Yc, self.section_Zc = extract_sections(sections)
        shear_center_ref = pd.read_excel(shear_center_ref)
        self.shear_center_y = (shear_center_ref['y'][0] * dimensions['b'] - self.dy)
        self.shear_center_z = (shear_center_ref['z'][0] * dimensions['h'] - self.dz)

        for node in self.nodes:
            node.y -= self.dy
            node.z -= self.dz

        self.force_dict['S_y'][1] = self.force_dict['S_y'][1] - self.dz
        self.force_dict['S_z'][1] = self.force_dict['S_z'][1] - self.dy
        self.qo = 0
        self.rate_of_twist_eqn = []


if __name__ == '__main__':
    shape = Shape()
    print(shape.dy, shape.dz, shape.Iy, shape.Iz, shape.Iyz)
    print(shape.nodes)
    print(shape.elements)
    print(shape.node_graph)
