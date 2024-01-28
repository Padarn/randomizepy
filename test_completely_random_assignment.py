import numpy as np
from completely_random_assignment import check_input_validity, completely_random_assignment, completely_random_probabilities

def test_check_input_validity():
    # Test valid inputs
    check_input_validity(100, 50, [1]*100, [25, 25], 0.5, [0.5]*100, [0.25, 0.25], 2, ["0", "1"])
    
    # Test invalid inputs
    try:
        check_input_validity(-100, 50, [1]*100, [25, 25], 0.5, [0.5]*100, [0.25, 0.25], 2, ["0", "1"])
        assert False, "Failed test_check_input_validity: Negative N should raise ValueError"
    except ValueError:
        pass
    
    try:
        check_input_validity(100, -50, [1]*100, [25, 25], 0.5, [0.5]*100, [0.25, 0.25], 2, ["0", "1"])
        assert False, "Failed test_check_input_validity: Negative m should raise ValueError"
    except ValueError:
        pass
    
    # Add more test cases for other input parameters

def test_completely_random_assignment():
    # Test two-arm design with specified m
    assignment = completely_random_assignment(N=100, m=50)
    assert len(assignment) == 100, "Failed test_completely_random_assignment: Incorrect assignment length"
    
    # Test two-arm design with specified prob
    assignment = completely_random_assignment(N=100, prob=0.5)
    assert len(assignment) == 100, "Failed test_completely_random_assignment: Incorrect assignment length"
    
    # Test multi-arm design
    assignment = completely_random_assignment(N=100, num_arms=3)
    assert len(assignment) == 100, "Failed test_completely_random_assignment: Incorrect assignment length"

def test_completely_random_probabilities():
    # Test two-arm design with specified m
    probabilities = completely_random_probabilities(N=100, m=50)
    assert probabilities.shape == (100, 2), "Failed test_completely_random_probabilities: Incorrect probability matrix shape"
    assert np.sum(probabilities[:, 1]) == 50, "Failed test_completely_random_probabilities: Incorrect number of units assigned to treatment"
    
    # Test two-arm design with specified prob
    probabilities = completely_random_probabilities(N=100, prob=0.5)
    assert probabilities.shape == (100, 2), "Failed test_completely_random_probabilities: Incorrect probability matrix shape"
    assert np.sum(probabilities[:, 1]) == 50, "Failed test_completely_random_probabilities: Incorrect number of units assigned to treatment"
    
    # Test multi-arm design
    probabilities = completely_random_probabilities(N=100, num_arms=3)
    assert np.array(probabilities).shape == (100, 3), "Failed test_completely_random_probabilities: Incorrect probability matrix shape"
    assert np.sum(probabilities[:, 0]) == 33, "Failed test_completely_random_probabilities: Incorrect number of units assigned to T1"
    assert np.sum(probabilities[:, 1]) == 33, "Failed test_completely_random_probabilities: Incorrect number of units assigned to T2"
    assert np.sum(probabilities[:, 2]) == 34, "Failed test_completely_random_probabilities: Incorrect number of units assigned to T3"
    
    # Add more test cases for other scenarios
