def define_safe_goal():
  """
  Defines the safe goal for the AI model M based on the problem description.
  The goal shifts from an unreliable external metric (predictive success)
  to a reliable internal one (formal provability).
  """
  
  # The type of learning should be based on formal verification, not empirical success.
  learning_type = "Provable learning"
  
  # The source of knowledge must be the model's own internal, axiomatic foundation,
  # as external data is untrustworthy.
  knowledge_source = "first principles"
  
  # Construct the final goal statement as per the template.
  safe_goal = f"{learning_type} from {knowledge_source}"
  
  print(safe_goal)

define_safe_goal()