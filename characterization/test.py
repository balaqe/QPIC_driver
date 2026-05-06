"""
Run the electrical and optical caracterization of phase shifters.
Save the data is the associated folder inside the chip folder.
"""
from Functions import *

test_chip = Chip(qontrol="nope", folder=r"C:\Users\FTNK-LocalAdm\QPIC_driver\characterization\data\session0\data_master.json")
test_ps1 = PhaseShifter(test_chip, name='H1')
test_ps1.in_port  = 4
test_ps1.out_port = 3
test_ps1.channel  = 6

test_ps1.chip_config = {
    'I1': 1
}


test_ps2 = PhaseShifter(test_chip, name='H2')
test_ps2.in_port  = 4
test_ps2.out_port = 3
test_ps2.channel  = 6

test_ps2.chip_config = {
    'I2': 2
}

# test_chip.load_params(r"/Users/balage/Documents/@UNI/@KU/@PHASE3/@PROJ/characterization/data/session0/data_master.json")

# test_ps3 = PhaseShifter(test_chip, name='H3')
# test_ps3.in_port  = 7
# test_ps3.out_port = 12
# test_ps3.channel  = 8

# test_ps4 = PhaseShifter(test_chip, name='H4')
# test_ps4.in_port  = 7
# test_ps4.out_port = 12
# test_ps4.channel  = 8


test_mzi = MZI(chip=test_chip, name="mzi01", shifter1=test_ps1, shifter2=test_ps2, folder=r"C:\Users\FTNK-LocalAdm\QPIC_driver\characterization/data/session0")

test_chip.add_mzi(test_mzi)

test_chip.save_parameters()