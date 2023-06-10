function generate_and_save_channel(numAntTx,numAntRx,numPath,numChn)
for ind = 1:numChn
    chn = generate_channel(numAntTx,numAntRx,numPath);
    save(['./data/channels_',num2str(numAntTx),'x',num2str(numAntRx),'/channel-',num2str(ind),'.mat'],'chn');
end
end