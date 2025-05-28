from .modeling_llama_megatron import (
    # original model with megatron
    ParallelLlamaModel,
    ParallelLlamaForCausalLM,
    # rmpad with megatron
    ParallelLlamaForCausalLMRmPad,
    ParallelLlamaForValueRmPad,
    # rmpad with megatron and pipeline parallelism
    ParallelLlamaForCausalLMRmPadPP,
    ParallelLlamaForValueRmPadPP)
