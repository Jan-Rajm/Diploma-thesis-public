timeStamp = string(datetime);
timeStamp = timeStamp.replace(":", "");
timeStamp = timeStamp.replace("-", "");
timeStamp = timeStamp.replace(" ", "");
m_simulationOutput = split(m_simulationOutput, ".");
m_simulationOutput = m_simulationOutput(1, 1);
for Path = inspected_paths
    m_simulationOutput = m_simulationOutput + Path.name;
end
save("stored_solutions_" + num2str(n_modes) + "modes_" + m_simulationOutput + "_" + timeStamp, "ALL_SOLUTIONS");