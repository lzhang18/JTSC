# Joint Time Series Chain (JTSC)

This is the source code for our SDM2022 paper: joint time series chain: detecting unusual evolving trend across time series. 

The main function is TSC_join. 

To run the demo example, please run Demo_JTSC after download the data. 

The data could be find in 
<a href="https://gmuedu-my.sharepoint.com/:u:/g/personal/lzhang18_gmu_edu/EQ6YXMdHxvxHvqqQTWr9WFgBKlb1PyCaIm8zCdv4c3G9eg?e=EC8dAZ"> here. </a>

The parameters are as follows: 

A: real-value array of normal data

B: real-value array of abnormal data 

f: file name for storing matrix profile

run_mp: 1 if we need to rerun matrix profile; 
        0 if we stored matrix profile before. In that case, we need to place a pre-computed mp in the same folder. 
        
SubseqLength: subsequence length of time series chain. Default value is set as 45.  

To run the demo in Matlab: 
>> Demo_JTSC

If you find the work or the data helpful, please cite the following: 

```bash
@inproceedings{zhang2022joint,
  title={Joint Time Series Chain: Detecting Unusual Evolving Trend across Time Series},
  author={Zhang, Li and Patel, Nital and Li, Xiuqi and Lin, Jessica},
  booktitle={Proceedings of the 2022 SIAM International Conference on Data Mining (SDM)},
  pages={208--216},
  year={2022},
  organization={SIAM}
}
```
