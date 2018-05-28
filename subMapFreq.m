function FreqData=subMapFreq(subCarrierData,subCars,FFtSize)
    a=reshape(subCarrierData,length(subCarrierData),1);
    FreqData=[0;a(subCars/2+1:end);zeros(FFtSize-subCars-1,1);a(1:subCars/2)];
    %return FreqData;
end