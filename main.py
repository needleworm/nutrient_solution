import system_simulator as ss

a = ss.Network("test_model.txt")
a.converge()
a.show_result()
