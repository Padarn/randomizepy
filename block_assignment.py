from typing import Optional, Union
import numpy as np


def block_random_assignment(
    blocks: np.ndarray, prob: Optional[float] = None, m: Optional[int] = None
) -> np.ndarray:
    """
    This function performs block random assignment.

    Args:
        blocks (np.ndarray): A numpy array that represents the block each unit belongs to.
        prob (Optional[float], optional): The probability of assigning a unit to treatment. Defaults to None.
        m (Optional[int], optional): The number of units to assign to treatment. Defaults to None.

    Returns:
        np.ndarray: A numpy array with the assignment of each unit. 1 represents treatment, 0 represents control.

    Note:
        If both `prob` and `m` are None, units are assigned to treatment or control randomly with equal probability.
        If `prob` is provided, approximately `prob` proportion of units in each block are assigned to treatment.
        If `m` is provided, exactly `m` units in each block are assigned to treatment.
    """
    unique_blocks = np.unique(blocks)
    assignment = np.empty_like(blocks)

    for block in unique_blocks:
        block_indices = np.where(blocks == block)[0]
        block_size = len(block_indices)

        if prob is not None:
            m = np.floor(block_size * prob)

        if m is not None:
            treatment_indices = np.random.choice(block_indices, int(m), replace=False)
            assignment[block_indices] = 0
            assignment[treatment_indices] = 1
        else:
            assignment[block_indices] = np.random.choice([0, 1], block_size)

    return assignment


def block_ra_probabilities(blocks=None, prob=None, prob_unit=None, prob_each=None, m=None, m_unit=None, block_m=None, block_m_each=None, block_prob=None, block_prob_each=None, num_arms=None, conditions=None, check_inputs=True):
  if check_inputs:
    .invoke_check(check_randomizr_arguments_new)
  else:
    N_per_block = pd.Series(blocks).value_counts().to_dict()
  
  block_spots = np.repeat(np.arange(len(blocks)), np.array([N_per_block[block] for block in blocks]))
  
  mapply_args = {
    'FUN': "complete_ra_probabilities",
    'N': N_per_block,
    'MoreArgs': {
      'conditions': conditions,
      'num_arms': num_arms,
      'check_inputs': False
    },
    'SIMPLIFY': False
  }
  
  prob_mat = block_ra_helper(blocks, prob, prob_unit, prob_each, m, m_unit, block_m, block_m_each, block_prob, block_prob_each, num_arms, N_per_block, mapply_args)
  
  prob_mat = pd.concat(prob_mat).sort_values(by=block_spots)
  
  return prob_mat


def block_ra_helper(blocks=None, prob=None, prob_unit=None, prob_each=None, m=None, m_unit=None, block_m=None, block_m_each=None, block_prob=None, block_prob_each=None, num_arms=None, N_per_block=None, mapply_args=None):
  if prob_unit is not None:
    block_prob = pd.Series(prob_unit).groupby(blocks).unique().to_dict()
  
  if m_unit is not None:
    block_m = pd.Series(m_unit).groupby(blocks).unique().to_dict()
  
  ret = {}
  
  if m is not None:
    ret['m'] = np.repeat(m, len(N_per_block))
  
  elif block_m is not None:
    ret['m'] = block_m
  
  elif block_prob is not None:
    ret['prob'] = block_prob
  
  elif prob_each is None and block_m_each is None:
    if prob is not None:
      prob_each = [1 - prob, prob]
    
    if prob_each is None:
      prob_each = [1 / num_arms] * num_arms
    
    ret['prob_each'] = [prob_each]
  
  elif block_m_each is not None:
    block_m_each_list = [block_m_each[i, :].tolist() for i in range(block_m_each.shape[0])]
    ret['m_each'] = block_m_each_list
  
  elif block_prob_each is not None:
    block_prob_each_list = [block_prob_each[i, :].tolist() for i in range(block_prob_each.shape[0])]
    ret['prob_each'] = block_prob_each_list
  
  return mapply(mapply_args, **ret)
