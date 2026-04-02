
## P0. `3,000,000 candidate sequences` / `500000 sequences based on low-frequency 5-mers` 这条链还有实质性问题，`random_02.R` 里存在硬伤


当前 `v1/xai/random_02.R` 里至少有两点需要修：

- 负样本部分单独生成了 `df2`，但写 `neg.csv` 时写出去的仍然是 `df`，不是 `df2`
- 前面虽然定义了 `bot_kmers`，但负样本生成部分当前仍然在用 `top_kmers`

这意味着当前仓库版本里，README 中下面这句还站不住：

- `500000 sequences based on low-frequency 5-mers`

同样，关于 `3,000,000 candidate sequences`，目前的问题也还没有完全闭合。  
如果你最终数据放在外部 `Dataset` 里，这当然可以，但至少要把下面几件事写清楚：

- `pos.csv` 的最终文件路径
- `neg.csv` 的最终文件路径
- 正样本到底是多少条
- 负样本到底是多少条
- 这两个文件各自的实际行数

这部分我建议不要再先写文稿数字、再回头解释仓库为什么暂时没对上。  
更稳妥的顺序应该是：

1. 先把脚本逻辑修对  
2. 再确认输出文件  
3. 最后再把 README 和文稿中的数字写死

## P0. `100,000 mutated sequences` 这一步当前仍然不是完整、独立、正式、可复跑脚本

你在 A2 里说：

- 输入文件是 `pos.csv`
- top10000 来自该文件
- 输出文件在 `Dataset/mut.csv`

但当前 `v1/xai/mut_04.R` 还没有达到这个状态。  
这里的问题已经不是“稍微不整洁”，而是脚本职责和文件路径都还不稳定。

当前版本至少存在下面这些问题：

- 读入路径写成了 `Dataset/pos.csv`
- `top10000` 是直接按行顺序取 `pos[1:10000, ]`，而不是明确从某个 ranked file 按某个排序规则选出
- 写出的文件名是 `mutant_sequence.csv`，不是 manifest 中写的 `../Dataset/mut.csv`
- 脚本后半段还混入了 mutant/raw 对比分布图分析，不是一个单一职责的 mutation generation script

这意味着 reviewer 仍然会追问：

- 你所谓的 `top10000` 到底是从哪个文件来的
- 是按 `score` 排序后的前 10000，还是只是当前文件前 10000 行
- mutation generation 的正式输出文件到底叫什么
- `100000` 是脚本真实产物，还是只是 `10000 x 10` 的理论配置值

## P0. ProteinX -> Rosetta -> final top5 这条链目前仍然是最需要澄清的一部分

这部分我还是建议你认真补，因为这里最容易被追问。

### 1. 书面答复和脚本参数目前还有直接冲突

你在 A3 里说：

- 每个序列保存 5 个构象

但 `v1/rosetta/run_rosetta_03.sh` 当前写的是：

- `NSTRUCT=20`

这两件事不是一定矛盾，但你需要明确解释它们分别对应哪一步。  
例如，如果意思是：

- ProteinX 每条序列先给 5 个模型
- Rosetta 再对每个模型做 20 个 decoys

那就把这件事写清楚。  
否则现在别人会很自然地看不懂“5 个构象”和 `20` 之间是什么关系。

### 2. `rosetta_05.R` 当前的 ranking 来源还不够清楚

现在这份脚本里，前面在算一套 ranking，后面又直接读取 `rank1.csv`，最后输出前 10 条。  
这里最容易引发的问题不是代码风格，而是结果来源会显得不够透明：

- 最终 ranking 到底来自脚本计算结果，还是来自外部 `rank1.csv`
- `result.csv` 的 top10 和文稿里的 final top5 到底是什么关系

如果你想表达的是：

- top10 用于展示
- top5 用于最终实验

这个逻辑没有问题。  
但最好在文件上也把这两件事分开，例如：

- 一个文件专门用于展示
- 一个文件专门用于最终筛选

否则现在这条链看起来还是有点混。

### 3. Rosetta manifest 现在还太简略

我看到你已经补了 `rosetta_06.json`，这方向是对的。  
但这份 manifest 现在还不足以回答 reviewer 最关心的问题：

- 输入序列文件是什么
- ProteinX 输出了什么
- Rosetta 实际吃进去的是什么
- 最终 top5 来自哪个 ranking 文件

这里不一定非要把所有中间结构文件都放进仓库，但至少应该把这些文件名和关系交代清楚。  
否则最终 top5 的 provenance 还是不够明白。

## P1. threshold 这次已经比上一版清楚很多，但还可以再收成一个更稳的来源链

最好让下游脚本直接从 manifest 读取 threshold，而不是在多个地方各自硬编码。


## P1. README 和 manifest 现在方向是对的，但还需要再对齐一次

这次 README / `v1/readme.md` 里的对照表是有帮助的，我建议保留。  
不过现在还差最后一步“对齐”：

- README 里写到的文件路径，要和仓库里实际存在的文件一致
- manifest 里写到的输出文件，要和外部 `Dataset` 中的文件名一致
- 如果数据放在 Google Drive，最好再多给一层说明：
  - 压缩包叫什么
  - 解压后目录叫什么
  - README / manifest 中的哪些路径对应到这个数据包中的哪些文件

这部分不是要你把所有东西都塞进仓库，而是要避免别人拿到仓库后还要自己猜路径关系。

## P1. `split_manifest_02.json` 可以再补一点关键信息

再补下面几个字段：

- `seed`
- `train_chr`
- `val_chr`
- `hold_chr`
