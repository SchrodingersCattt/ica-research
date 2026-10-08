import json
import glob
import os

pretraining_trn = glob.glob("../../EMs-v1/train_cleaned/*/")
finetuning_trn = glob.glob("new_data/*/*/trn_*/")
trn = pretraining_trn + finetuning_trn
val = glob.glob("new_data/*/*/val_*/")
pretrain_idx = len(pretraining_trn)
finetune_idx = len(finetuning_trn)

def get_sys_prob(pretrain_idx, finetuning_trn):
    sys_prob = []
    for i in range(pretrain_idx):
        sys_prob.append(0.5/pretrain_idx)
    for i in range(finetune_idx-1):
        sys_prob.append(0.5/finetune_idx)
    sys_prob += [1 - sum(sys_prob)]  ## force the overall prob to be 1
    return sys_prob

print(len(trn), pretrain_idx, finetune_idx)
config = "template.json"

with open(config) as f:
    config_dict = json.load(f)
config_dict["training"]["training_data"]["systems"] = [os.path.abspath(x) for x in trn]
config_dict["training"]["validation_data"]["systems"] = [os.path.abspath(x) for x in val]
sys_prob = get_sys_prob(pretrain_idx, finetuning_trn)

print(sys_prob, sum(sys_prob))
config_dict["training"]["training_data"]["sys_probs"] = sys_prob
with open(config, "w") as f:
    json.dump(config_dict, f, indent=4)
