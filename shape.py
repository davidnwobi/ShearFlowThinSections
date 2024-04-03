from data_setup import extract_nodes_and_elements, convert_nodes_and_elements_to_shape, calculate_beam_properties, \
    extract_forces, extract_sections


class Shape:
    def __init__(self, nodes="Nodes.xlsx", elements="Elements.xlsx", forces="Forces.xlsx"):
        nodes_df, elements_df, self.nodes, self.elements, self.node_graph = extract_nodes_and_elements(nodes,
                                                                                                       elements)
        shape = convert_nodes_and_elements_to_shape(nodes_df, elements_df)
        self.force_dict = extract_forces(forces)
        self.dy, self.dz, self.Ac, self.Iy, self.Iz, self.Iyz = calculate_beam_properties(shape)
        self.ey, self.ez = 0, 0
        self.sections, self.areas = extract_sections()

        for node in self.nodes:
            node.y -= self.dy
            node.z -= self.dz

        self.force_dict['S_y'][1] = self.force_dict['S_y'][1] - self.dz
        self.force_dict['S_z'][1] = self.force_dict['S_z'][1] - self.dy
        self.qo = 0


if __name__ == '__main__':
    shape = Shape()
    print(shape.dy, shape.dz, shape.Iy, shape.Iz, shape.Iyz)
    print(shape.nodes)
    print(shape.elements)
    print(shape.node_graph)
