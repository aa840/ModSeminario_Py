def unit_vector_N(u_BC, u_AB):
    # Calculates unit normal vector which is perpendicular to plane ABC
    import numpy as np

    cross_product = np.cross(u_BC, u_AB)
    norm_u_N = np.linalg.norm(cross_product)
    u_N = cross_product / norm_u_N
    return u_N
