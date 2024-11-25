function runParallelJob(inputFilePath, taskIndex, jobNumber, parameterName, parameterValues)

newFolderPathList = generateDataFolders(inputFilePath, taskIndex, jobNumber, parameterName, eval(parameterValues));

while (~exist([newFolderPathList{taskIndex} 'inputFile_' num2str(jobNumber) '.json'], 'file'))
end

main(newFolderPathList{taskIndex}, ['inputFile_' num2str(jobNumber) '.json'], 0, 1)

end
