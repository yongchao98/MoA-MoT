import os
import csv
import json
import reasoning_gym
import numpy as np

seed = 32
task_size = 1000

available_datasets = [
    # 'ab',
    # 'acre',
    # 'advanced_geometry',
    # 'aiw',
    # 'arc_1d',
    # 'arc_agi',
    # 'base_conversion',
    # 'basic_arithmetic',
    # 'bf',
    # 'binary_alternation',
    # 'binary_matrix',
    # 'bitwise_arithmetic',
    # 'caesar_cipher',
    # 'calendar_arithmetic',
    # 'chain_sum',
    'circuit_logic',
]

output_dir = './dataset_gather/reasoning_gym'
os.makedirs(output_dir, exist_ok=True)

def convert_to_serializable(obj):
    if isinstance(obj, (np.integer,)):
        return int(obj)
    elif isinstance(obj, (np.floating,)):
        return float(obj)
    elif isinstance(obj, (np.ndarray,)):
        return obj.tolist()
    elif hasattr(obj, 'item'):  # for torch.Tensor scalar
        return obj.item()
    else:
        return str(obj)  # fallback for unknown types

for dataset in available_datasets:
    data = reasoning_gym.create_dataset(dataset, size=task_size, seed=seed)
    output_path = os.path.join(output_dir, f'{dataset}.csv')

    with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['ID', 'dataset', 'question', 'answer', 'full_data'])
        writer.writeheader()

        for i, x in enumerate(data):
            # print(f'x: {x}')
            # Optional: validate that the answer is correct
            # assert data.score_answer(answer=x['answer'], entry=x) == 1.0

            question = x['question'] + '\nOutput final answer with the format <<<answer>>>'
            if dataset == 'arc_agi':
                question = x['question'].replace('Your final answer should just be the text output grid itself.',
                                                 'Your final answer should be the output grid enclosed in triple angle brackets, like this: <<<output grid>>>')
            elif dataset == 'bf':
                question = x['question'].replace('Respond only with the exact output of the program.',
                                                 'Respond only with the exact output of the program enclosed in triple angle brackets, like this: <<<output>>>.')
            elif dataset == 'boxnet':
                question = x['question'] + '\nOutput action plan enclosed in triple angle brackets, like this <<<action plan>>>'

            writer.writerow({
                'ID': i,
                'dataset': dataset,
                'question': question,
                'answer': x['answer'],
                'full_data': json.dumps(x, indent=4, ensure_ascii=False, default=convert_to_serializable)  # serialize full x dictionary as string
            })

    print(f"Saved {task_size} tasks from '{dataset}' to '{output_path}'")
# print(reasoning_gym.factory.DATASETS)