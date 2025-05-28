import os
from dataclasses import dataclass
from verl.single_controller.base.worker import Worker, DistRankInfo, DistGlobalInfo


class MegatronWorker(Worker):

    def __init__(self, cuda_visible_devices=None) -> None:
        super().__init__(cuda_visible_devices)

    def get_megatron_global_info(self):
        from megatron.core import parallel_state as mpu
        tp_size = mpu.get_tensor_model_parallel_world_size()
        dp_size = mpu.get_data_parallel_world_size()
        pp_size = mpu.get_pipeline_model_parallel_world_size()
        info = DistGlobalInfo(tp_size=tp_size, dp_size=dp_size, pp_size=pp_size)
        return info

    def get_megatron_rank_info(self):
        from megatron.core import parallel_state as mpu
        tp_rank = mpu.get_tensor_model_parallel_rank()
        dp_rank = mpu.get_data_parallel_rank()
        pp_rank = mpu.get_pipeline_model_parallel_rank()
        info = DistRankInfo(tp_rank=tp_rank, dp_rank=dp_rank, pp_rank=pp_rank)
        return info