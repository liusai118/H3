## P0. 负样本生成逻辑仍然和文稿/README 的说法不一致，这不是表述问题，而是 decoy 集合定义本身没有对上


- `df2` 现在已经被写到 `neg.csv`，这点确实修复了，
但更关键的问题还在：当前 `v1/xai/random_02.R` 里，负样本生成部分仍然在用 `top_kmers`，而不是前面定义的 `bot_kmers`。

这不是 README 文案精不精确的问题，而是**负样本 / decoy 集合到底是按什么规则构造的**没有对上。  
一旦 decoy 的生成逻辑和 Methods 中写的不一致，reviewer 会直接质疑：

- 你后续关于正负候选分离的叙述，到底基于哪一套 decoy 定义
- 这个 negative set 是不是你文稿里声称的 low-frequency 5-mer decoy
- 如果不是，那后面的比较和解释是否仍然成立

当前证据非常具体：

- `v1/xai/random_02.R` 里定义了 `bot_kmers <- kmers[925:1024]`
- 但生成 `df2` 时传入的仍然是 `top_kmers`
- README 仍然写的是 low-frequency 5-mers

所以这条我建议直接按下面的标准修：

1. 先把 `neg.csv` 的生成逻辑改成真正使用 `bot_kmers`  
2. 再确认 `neg.csv` 的实际行数  
3. 最后再更新 README 和文稿中的数字/描述  


## P0. `100,000 mutated sequences` 这一步当前仍然不能作为独立可复跑证据，因为脚本本身还没有锁成一个稳定的正式产物生成步骤


目前还存在三个具体问题：

### 1. 当前脚本按字面形式不能直接运行

`v1/xai/mut_04.R` 第一行现在仍然是：

```r
pos <- read.csv(Dataset/pos.csv)
```
这行代码按当前形式就不能直接执行。  
`read.csv()` 这里需要字符串路径；现在没有引号，脚本本身并不能按当前提交版本独立运行。

### 2. 输出文件名仍然没有统一

当前几处说法仍然不一致：

- `v1/xai/mut_04.R` 写出的是 `mut.csv`
- `v1/script_description/mutation_manifest_05.json` 写的是 `../Dataset/mut.csv`
- `README.md` 写的是 `./Dataset/mutant.csv`

这不是小问题。  
因为 reviewer 在这里要找的是：**文稿中的 100,000 mutated sequences 到底对应哪个最终文件**。  
只要文件名还在三处各不相同，这条证据链就还没锁住。

### 3. mutation generation 和后续分析仍然混在同一个脚本里

当前 `mut_04.R` 里同时在做：

- mutant 生成
- mutant 评分后的合并分析
- `ProteinX.csv` 导出

这会让 reviewer 很难判断：

- 哪一步才是正式的 mutation generation
- 哪个文件才是这一步的主输出
- `ProteinX.csv` 到底是从 mutation generation 直接得到，还是从后续分析再次筛选得到

如果你论文里把 `100,000 mutated sequences` 和后面的 ProteinX / Rosetta 链条当成连续证据，那么这里的正式输出边界必须更清楚。

所以这条我建议按 reviewer 能接受的最低标准处理：

1. 先把 `mut_04.R` 变成一个只负责生成 `mut.csv` 的脚本  
2. 明确输入文件、排序字段、top10000 选择规则、输出文件名  
3. 把后续 analysis / ProteinX 导出放到单独脚本  

不然 reviewer 会自然怀疑：  
你这 10 万条 mutant 到底是已经真实生成并固定下来的结果，还是脚本里若干步混在一起后的中间状态。

## P0. final top5 的 provenance 现在仍然没有完全闭合；这是我认为当前最接近“会被 reviewer 直接追住不放”的问题

这条我仍然建议你优先处理。  
因为 ProteinX -> Rosetta -> final top5 如果是论文主筛选链条的一部分，那么 reviewer 真正要看的不是你有没有做图，而是：

- 最终实验采用的 top5 到底来自哪个具体文件
- 这个文件又是从哪个 ranking 步骤产生的
- 中间有没有手工替换或后验挑选

### 1. 用于展示和用于实验的文件已经分开，但 README / manifest 还没有跟上

- `README.md` 仍然写的是 `result.csv`
- `v1/script_description/rosetta_vis_07.json` 也仍然写的是 `result.csv`

这意味着 reviewer 现在仍然看不清：

- 文稿里的 final top5 对应的是 `exp.csv`，还是别的文件
- README 指向的 `result.csv` 和现在代码里的 `exp.csv` / `exh.csv` 到底哪个才是最终版本

只要 final top5 不能唯一映射到一个具体文件，这条链就是不闭合的。

### 2. Rosetta manifest 仍然太笼统，不能回答 final top5 的来源问题

当前 `v1/script_description/rosetta_run_06.json` 里写的是：

- `"input_file": "output of ProteinX"`

这对审稿是远远不够的。  
因为 reviewer 关心的不是“原则上 Rosetta 的输入来自 ProteinX”，而是：

- 具体是哪一个文件
- 一共多少条序列
- 每条序列保留了几个 ProteinX 模型
- Rosetta 吃进去的是哪个目录/哪些结构文件
- top5 最终来自哪个 ranking 文件

如果这些问题回答不具体，reviewer 很容易把这条链看成：

- 中间做了很多结构预测
- 但最终表里的 5 条序列来源并不完全透明

### 3. 这条链里还有一个参数解释需要锁死

你现在的文字解释是：

- ProteinX 每条序列保留 5 个构象
- Rosetta 再对每个构象分别优化 20 次

这逻辑是可以成立的。  
但请一定把这件事在 manifest 或 release 说明里写成明确流程，而不要只保留在答复文本里。  
因为只要 reviewer 只能在问答里看到这层关系，而代码/manifest 里没有对应说明，他们还是会认为 provenance 不够稳定。

所以这条我建议你最终至少做到：

1. README 只保留一个 final top5 文件名  
2. manifest 里明确写出 ProteinX 输入、ProteinX 输出、Rosetta 输入、Rosetta score 文件、final top5 文件  
3. 文稿中的 final top5 明确指向这个文件

否则 reviewer 会很自然地追问：**你文章里的 final five 到底是代码固定产物，还是后期整理出来的表。**

## P0. 当前 release 对“10 个 replicate”这件事的支持是不完整的；`rep1` 数据缺失会直接伤到你 10-rep 结果的可复现性

当前仓库里`rep1` 对应的数据文件并不存在。

你现在不只是“少放了一个例子文件”，而是：

- 你在文稿、threshold manifest、run notebook 里都在使用 10-rep 的证据
- 但当前 release 并不能把这 10-rep 的来源完整地交付出来

reviewer 一旦看到这点，会直接问：

- 10 个 rep 到底是不是真的都在当前 release 中可复查
- `rep1` 对应的 threshold 和训练结果，当前是否还能独立重现
- 如果 `rep1` 数据不在 release 中，为什么你还把它算在 10-rep 汇总里

这条我建议不要轻描淡写。  
因为这不是“路径还可以再整理”，而是**10-replicate evidence 在当前 release 中缺了一个实际输入分支**。

在这一步补齐之前，任何涉及：

- 10-rep threshold 汇总
- 10-rep model performance
- 10-rep candidate scoring

的叙述，都会被 reviewer 认为当前 archive 支撑不完整。

## P1. 这条 pipeline 目前仍然有明显的 self-reinforcement 风险；如果文稿继续把它写成独立 discovery 证据，reviewer 很可能会卡

- 你先用同一模型家族在 dev 序列上做打分 / IG / kmer 选择
- 再用这些高分 kmer 去拼接生成候选序列
- 然后再用同一模型家族去对这些候选序列重新打分和排序

当前代码链条是很清楚的：

- `v1/xai/kmer_choose_01.ipynb` 从 model-scored top sequences 中统计高分 k-mer
- `v1/xai/random_02.R` 用这些 top k-mers 生成 candidate sequences
- `v1/xai/score_03.py` 再用同一批模型对 `pos.csv` 重新打分

从 reviewer 角度看，这会带来一个非常自然的质疑：

- 这些 top candidates 到底是“被独立证据支持”的候选，还是“被同一模型偏好强化出来”的候选

这条不一定单独导致拒稿；  
但如果文稿把这部分写成较强的 discovery 结论，而不是更谨慎地写成：

- `model-guided prioritization`

那么 reviewer 很可能会认为这里存在 circularity / self-reinforcement 风险。

所以这条我建议至少做下面之一：

1. 在 Methods / Discussion 里明确承认这是 model-guided candidate generation + rescoring 的闭环流程  
2. 增加一个更正交的 ranking / filtering / validation 证据，来说明最终 top hits 不是单纯由同一模型偏好反复放大的结果  

如果这部分继续写得过满，reviewer 很可能不会接受“这些 top hits 代表独立发现”。
