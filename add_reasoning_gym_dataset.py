import os
import csv
import json
import reasoning_gym

seed = 32
task_size = 1000

available_datasets = [
    'ab'
]

output_dir = './dataset_gather/reasoning_gym'
os.makedirs(output_dir, exist_ok=True)

for dataset in available_datasets:
    data = reasoning_gym.create_dataset(dataset, size=task_size, seed=seed)
    output_path = os.path.join(output_dir, f'{dataset}.csv')

    with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=['ID', 'dataset', 'question', 'answer', 'full_data'])
        writer.writeheader()

        for i, x in enumerate(data):
            # Optional: validate that the answer is correct
            # assert data.score_answer(answer=x['answer'], entry=x) == 1.0

            writer.writerow({
                'ID': i,
                'dataset': dataset,
                'question': x['question'] + '\nOutput final answer with the format <<<answer>>>',
                'answer': x['answer'],
                'full_data': json.dumps(x, indent=4, ensure_ascii=False)  # serialize full x dictionary as string
            })

    print(f"Saved {task_size} tasks from '{dataset}' to '{output_path}'")