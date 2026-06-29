import numpy as np
import poreana as pa


def test_utils():
    pa.utils.mkdirp("output/test")

    # .obj round-trip
    file_link = "output/test/test.obj"
    pa.utils.save([1, 1, 1], file_link)
    assert pa.utils.load("output/test/test.txt", file_type="DOTA") is None
    assert pa.utils.load(file_link) == [1, 1, 1]

    # .yml round-trip (also exercises _to_python conversion)
    yml_link = "output/test/test.yml"
    pa.utils.save({"a": 1, "b": [2.0, 3.0]}, yml_link)
    loaded = pa.utils.load(yml_link)
    assert loaded["a"] == 1
    assert loaded["b"] == [2.0, 3.0]

    # unit-conversion helpers
    assert round(pa.utils.mumol_m2_to_mols(3, 100), 4) == 180.66
    assert round(pa.utils.mols_to_mumol_m2(180, 100), 4) == 2.989
    assert round(pa.utils.mmol_g_to_mumol_m2(0.072, 512), 2) == 0.14
    assert round(pa.utils.mmol_l_to_mols(30, 1000), 4) == 18.066
    assert round(pa.utils.mols_to_mmol_l(18, 1000), 4) == 29.8904

    pa.utils.toc(pa.utils.tic(), message="Test", is_print=True)
    assert round(pa.utils.toc(pa.utils.tic(), is_print=True)) == 0


def test_geometry():
    vec_a = [1, 1, 2]
    vec_b = [0, 3, 2]
    assert round(pa.geom.dot_product(vec_a, vec_b), 4) == 7
    assert round(pa.geom.length(vec_a), 4) == 2.4495
    assert [round(x, 4) for x in pa.geom.vector(vec_a, vec_b)] == [-1, 2, 0]
    assert pa.geom.vector([0, 1], [0, 0, 0]) is None
    assert [round(x, 4) for x in pa.geom.unit(vec_a)] == [0.4082, 0.4082, 0.8165]
    assert [round(x, 4) for x in pa.geom.cross_product(vec_a, vec_b)] == [-4, -2, 3]
    assert round(pa.geom.angle(vec_a, vec_b), 4) == 37.5714
