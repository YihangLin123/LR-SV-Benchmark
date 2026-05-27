# LR-SV-Benchmark

## Contents
1. [Introduction](#introduction)
2. [Datasets](./Docs/Datasets.md)
3. [SV Eallers](./Docs/SVCallers.md)
4. [SV Evaluators](./Docs/SVEvaluators.md)
5. [Detailed running command]
6. [Data availability](#data-availability)

## Introduction

In this study, we developed a unified evaluation framework encompassing both germline and somatic structural variant (SV) detection to systematically assess the downstream performance of various callers. 

A total of 14 SV detection tools were benchmarked and stratified into three functional categories based on their application scenarios:
* Germline-specific tools (9 callers)
* Somatic-dedicated algorithms (2 callers)
* Versatile callers (3 callers, compatible with both modes)

This comprehensive framework provides reliable recommendations for tool selection across diverse research contexts while promising to guide future algorithmic development. The overall evaluation workflow of this study is illustrated in the figure below:

<p align="center">
  <img src="./Docs/Fig%201.png" alt="Overall Evaluation Workflow" width="80%">
</p>

## Data availability
