def check_randomizr_arguments(N=None, prob=None, prob_unit=None, m=None, m_unit=None, m_each=None, prob_each=None, blocks=None, block_m=None, block_m_each=None, block_prob=None, block_prob_each=None, clusters=None, num_arms=None, simple=None, conditions=None, **kwargs):
  if clusters is not None:
    N = len(set(clusters))
  
  if blocks is not None:
    if simple:
      raise ValueError("You can't specify `simple` when using blocked assignment.")
    
    if clusters is not None:
      N_per_block = {b: len(set(c)) for b, c in zip(blocks, clusters)}
      if any(sum(table(blocks, clusters) > 0) > 1):
        raise ValueError("All units within a cluster must be in the same block.")
    else:
      N_per_block = {b: blocks.count(b) for b in blocks}
    
    if N is None:
      N = sum(N_per_block.values())
    elif N != sum(N_per_block.values()):
      raise ValueError("N should equal the length of blocks.")
    
    N_blocks = len(N_per_block)
    blocks = {"N_per_block": N_per_block, "N_blocks": N_blocks}
  
  if N is None:
    raise ValueError("N, blocks or clusters must be specified.")
  
  if not isinstance(N, int) or N <= 0:
    raise ValueError("N must be a positive integer scalar.")
  
  if len(set(conditions)) != len(conditions):
    raise ValueError("You must supply unique values to conditions.")
  
  if num_arms is not None and conditions is not None and num_arms != len(conditions):
    raise ValueError("If both conditions and num_arms are specified, the length of conditions must be equal to num_arms.")
  
  conflict_args = ["prob", "prob_each", "prob_unit", "m", "m_unit", "m_each", "block_prob", "block_prob_each", "block_m", "block_m_each"]
  specified_args = {arg: value for arg, value in kwargs.items() if arg in conflict_args}
  
  if len(specified_args) > 1:
    raise ValueError("Please specify only one of " + " and ".join(specified_args.keys()) + ".")
  elif len(specified_args) == 1:
    arg = next(iter(specified_args))
    arg_block = arg.startswith("block_")
    arg_each = arg.endswith("_each")
    
    if arg_block and blocks is None:
      raise ValueError("Specified `" + arg + "` but blocks is None.")
    
    if simple and "prob" not in arg:
      raise ValueError("You can't specify `" + arg + "` when simple = True.")
    
    if num_arms is not None and conditions is not None:
      _check_ra_arg_num_arms_conditions(arg, arg_block, arg_each, specified_args[arg], num_arms, conditions)
    
    _check_ra[arg](N, blocks, clusters, num_arms, conditions, simple, specified_args[arg])
  
  if num_arms is None:
    if conditions is not None:
      num_arms = len(conditions)
    elif len(specified_args) == 0:
      num_arms = 2
    elif not arg_each:
      num_arms = 2
    elif not arg_block and not isinstance(specified_args[arg], list):
      num_arms = len(specified_args[arg])
    else:
      num_arms = len(specified_args[arg])
    
    if num_arms == 2 and conditions is None:
      conditions = [0, 1]
  
  if conditions is None:
    conditions = ["T" + str(i) for i in range(1, num_arms + 1)]
  
  ret = {
    "num_arms": num_arms,
    "conditions": conditions,
    "N_per_block": blocks.get("N_per_block")
  }
  
  return ret


class CheckRA:
    def __init__(self):
        pass

    def check_ra_arg_num_arms_conditions(arg, arg_block, arg_each, value, num_arms, conditions):
        if not arg_each:
          w = 2
          num_arms_fmt = "If {} and num_arms are both specified, {} must be length 2 and num_arms must be 2."
          conditions_fmt = "If {} and conditions are both specified, {} and conditions must both be length 2."
        elif not arg_block:
          w = len(value)
          num_arms_fmt = "If {} and num_arms are both specified, the length of {} must be equal to num_arms."
          conditions_fmt = "If {} and conditions are both specified, the length of {} must be equal to the length of conditions."
        else:
          w = value.shape[1]
          num_arms_fmt = "If {} and num_arms are both specified, the number of columns of {} must be equal to num_arms."
          conditions_fmt = "If {} and conditions are both specified, the length of conditions must be equal to the number of columns of {}."
        
        if num_arms is not None and num_arms != w:
          raise ValueError(num_arms_fmt.format(arg, arg))
        
        if conditions is not None and len(conditions) != w:
          raise ValueError(conditions_fmt.format(arg, arg))

    def prob(N, blocks, clusters, num_arms, conditions, simple, prob):
        if any(p > 1 or p < 0 for p in prob):
          raise ValueError("The probability of assignment to treatment must be between 0 and 1.")
        if len(prob) not in [1]:
          raise ValueError("`prob` must be of length 1.")

    def prob_unit(self, N, blocks, clusters, num_arms, conditions, simple, prob_unit):
        if any(p > 1 or p < 0 for p in prob_unit):
          raise ValueError("The probability of assignment to treatment must be between 0 and 1.")
        
        if blocks is None:
          if not simple:
            if not all(p == prob_unit[0] for p in prob_unit):
              raise ValueError("In a complete random assignment design, `prob_unit` must be the same for all units.")
        else:
          if not all(p == prob_unit[0] for p in prob_unit):
            raise ValueError("In a block random assignment design, `prob_unit` must be the same for all units within the same block.")
        
        if clusters is not None:
          if len(prob_unit) != len(clusters):
            raise ValueError("`prob_unit` must be of length N.")
        else:
          if len(prob_unit) != N:
            raise ValueError("`prob_unit` must be of length N.")
          
    def prob_each(self, N, blocks, clusters, num_arms, conditions, simple, prob_each):
        if any(p > 1 or p < 0 for p in prob_each):
          raise ValueError("The probabilities of assignment to any condition may not be greater than 1 or less than zero.")
        if isinstance(prob_each, list):
          if sum(prob_each) != 1:
            raise ValueError("The sum of the probabilities of assignment to each condition (prob_each) must equal 1 for each obs.")
        elif isinstance(prob_each, np.ndarray):
          if any(np.sum(prob_each, axis=1) != 1):
            raise ValueError("The sum of the probabilities of assignment to each condition (prob_each) must equal 1 for each obs.")
        function(N,
                blocks,
                clusters,
                num_arms,
                conditions,
                simple,
                prob_each) {
          if (any(prob_each > 1 | prob_each < 0)) {
            stop(
              "The probabilties of assignment to any condition may not be greater than 1 or less than zero.",
              call. = FALSE
            )
          }
          if (is.vector(prob_each)) {
            if (all.equal(sum(prob_each), 1) != TRUE) {
              stop(
                "The sum of the probabilities of assignment to each condition (prob_each) must equal 1 for each obs.",
                call. = FALSE
              )
            }
          }
          else if (is.matrix(prob_each)) {
            if (any(sapply(rowSums(prob_each), function(x)
              all.equal(x, 1) != TRUE))) {
              stop(
                "The sum of the probabilities of assignment to each condition (prob_each) must equal 1 for each obs.",
                call. = FALSE
              )
            }
            if (!nrow(prob_each) %in% c(1, N)) {
              stop("`prob_each` must have either 1 or N rows.", call. = FALSE)
            }
          }
          else
            stop("`prob_each` must be a vector or matrix.", call. = FALSE)
          
            
  

    def m(self, N, blocks, clusters, num_arms, conditions, simple, m):
        if any(p > 1 or p < 0 for p in prob_each):
          raise ValueError("The probabilities of assignment to any condition may not be greater than 1 or less than zero.")
        if isinstance(prob_each, list):
          if sum(prob_each) != 1:
            raise ValueError("The sum of the probabilities of assignment to each condition (prob_each) must equal 1 for each obs.")
        elif isinstance(prob_each, np.ndarray):
          if any(np.sum(prob_each, axis=1) != 1):
            raise ValueError("The sum of the probabilities of assignment to each condition (prob_each) must equal 1 for each obs.")

    def m_unit(self, N, blocks, clusters, num_arms, conditions, simple, m_unit):
         # if it's complete
        if blocks is None:
          if not all(x == m_unit[0] for x in m_unit):
            raise ValueError("In a complete random assignment design, `m_unit` must be the same for all units.")
        else:
          # if it's blocked
          for block in blocks:
            if not all(x == m_unit[block[0]] for x in m_unit[block]):
              raise ValueError("In a block random assignment design, `m_unit` must be the same for all units within the same block.")
        
        # if it's clustered:
        if clusters is not None:
          if len(m_unit) != len(clusters):
            raise ValueError("`m_unit` must be of length N.")
        else:
          if len(m_unit) != N:
            raise ValueError("`m_unit` must be of length N.")

    def m_each(self, N, blocks, clusters, num_arms, conditions, simple, m_each):
        # similar logic as in your R function

    def block_prob(self, N, blocks, clusters, num_arms, conditions, simple, block_prob):
        # similar logic as in your R function

    def block_prob_each(self, N, blocks, clusters, num_arms, conditions, simple, block_prob_each):