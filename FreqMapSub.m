function subData=FreqMapSub(FreqData,subCars)
    a=reshape(FreqData,length(FreqData),1);
    subData=[a(end-(subCars/2)+1:end);a(2:subCars/2+1)];
    %return subData;
end