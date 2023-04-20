function generate_and_save_channel_mu(numUser,numAntTx,numAntRx,numPath,numChn,type)
% GENERATE_AND_SAVE_CHANNEL_MU generates channels for multiple users and
% save it in a predefined place.
%
% Inputs:
% - numUser [scalar]: the number of users
% - numAntTx [scalar]: the number of transmit antennas
% - numAntRx [scalar]: the number of receiver antennas per user
% - numPath [scalar]: the number of paths in the channel
% - numChn [scalar]: the number of channels to be generated
% - type [string]: the type of function output, 'combined' or 'separated'.
%     The default type is 'combined'. The combined channel for all K users
%     is
%                       H = [h1,h2,...,hK]^H,
%     where hk is the channel
%
% Outputs:
% 


for ind = 1:numChn

    % Generate Multiuser Channel
    if type == "combined"
        chn = zeros(numAntRx*numUser, numAntTx);
        for iUser = 1:numUser
            chn((iUser-1)*numAntRx+1:iUser*numAntRx, :) = generate_channel(numAntTx,numAntRx,numPath);
        end
    elseif type == "separated"
        chn = zeros(numAntRx,numAntTx,numUser);
        for iUser = 1:numUser
            chn(:,:,iUser) = generate_channel(numAntTx,numAntRx,numPath);
        end
    else
        error("Not supported type of output variable!");
    end
    
    save(['./data/channels_mu_',num2str(numAntTx),'x',num2str(numAntRx),'/channel-',num2str(ind),'.mat'],'chn');
end
end