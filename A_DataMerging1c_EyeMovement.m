%% A_DataMerging1c
% 2018.10.8

search_result(:,1)=k+cond2-1;
mem_ans=[sub_answer(:,1) ref_answer(:,1)];
memory_dev=2*abs(mem_ans(:,1)-mem_ans(:,2));
% Calculate deviation
for k=1:length(memory_dev)
    if memory_dev(k)>180
        memory_dev(k)=360-memory_dev(k);
    end
end
memory_dev(memory_dev>180)=360-memory_dev(memory_dev>180);
seq_dis_pos_all=repmat(seq_dis_pos',8,1);
seq_con=seq_con';
seq_gap=seq_gap';
seq_tar_pos=seq_tar_pos';
test_ans=test_ans';
search_m=table(search_result, seq_con, seq_gap, seq_tar_pos,...
    seq_dis_pos_all, mem_ans, memory_dev, test_ans);