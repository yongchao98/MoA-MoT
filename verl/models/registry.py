import importlib
from typing import List, Optional, Type

import torch.nn as nn

# Supported models using HF Rmpad
# TODO(sgm): HF may supported more than listed here, we should add more after testing
from transformers import LlamaConfig, MistralConfig, GemmaConfig, Qwen2Config

_REOVEPAD_MODELS = {'llama': LlamaConfig, 'mistral': MistralConfig, 'gemma': GemmaConfig, 'qwen2': Qwen2Config}


def check_model_support_rmpad(model_type: str):
    assert isinstance(model_type, str)
    if not model_type in _REOVEPAD_MODELS.keys():
        raise ValueError(f"Model architecture {model_type} is not supported for now. "
                         f"RMPad supported architectures: {_REOVEPAD_MODELS.keys()}."
                         f"Please set `use_remove_padding=False` in the model config.")


# Supported models in Megatron-LM
# Architecture -> (module, class).
_MODELS = {
    "LlamaForCausalLM":
        ("llama", ("ParallelLlamaForCausalLMRmPadPP", "ParallelLlamaForValueRmPadPP", "ParallelLlamaForCausalLMRmPad")),
    "MistralForCausalLM": ("mistral", ("ParallelMistralForCausalLMRmPadPP", "ParallelMistralForValueRmPadPP",
                                       "ParallelMistralForCausalLMRmPad"))
}


# return model class
class ModelRegistry:

    @staticmethod
    def load_model_cls(model_arch: str, value=False) -> Optional[Type[nn.Module]]:
        if model_arch not in _MODELS:
            return None

        megatron = "megatron"

        module_name, model_cls_name = _MODELS[model_arch]
        if not value:  # actor/ref
            model_cls_name = model_cls_name[0]
        elif value:  # critic/rm
            model_cls_name = model_cls_name[1]

        module = importlib.import_module(f"verl.models.{module_name}.{megatron}.modeling_{module_name}_megatron")
        return getattr(module, model_cls_name, None)

    @staticmethod
    def get_supported_archs() -> List[str]:
        return list(_MODELS.keys())
