def get_weight_loader(arch: str):
    from verl.models.llama.megatron.checkpoint_utils.llama_loader import load_state_dict_to_megatron_llama
    _MODEL_WEIGHT_MEGATRON_LOADER_REGISTRY = {'LlamaForCausalLM': load_state_dict_to_megatron_llama}

    if arch in _MODEL_WEIGHT_MEGATRON_LOADER_REGISTRY:
        return _MODEL_WEIGHT_MEGATRON_LOADER_REGISTRY[arch]
    raise ValueError(f"Model architectures {arch} are not supported for now. "
                     f"Supported architectures: {_MODEL_WEIGHT_MEGATRON_LOADER_REGISTRY.keys()}")
