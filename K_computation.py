import simulation

om_1_set = []
B_set = [1.1, 1.3, 1.5, 1.7]

for B in B_set:
    s = simulation.Simulation()
    s.B = B
    s.T = 2**14
    s.solve()
    s.transform(R=True)
    s.find_peaks()
    om_1_set.append(s.peaks[0])

om_b = om_1_set[0]

K_set = [om_b/om for om in om_1_set]

print(K_set)
