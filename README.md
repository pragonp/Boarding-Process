# Airbus A320 Boarding Simulation
Welcome to my first coding project from 2018! This repository contains the original code and serves as a snapshot of my early programming work, focused on simulating the boarding process for an Airbus A320.

## Overview
This project is based on the Social Force Model developed by Dirk Helbing and Péter Molnar. It was suggested to me by Dr. Matthew Borg, and I was particularly intrigued by the intersection of technology and social behaviors. The study aims to model agent interactions during boarding and evaluate different strategies.

This project was recognized by the University Head of Department at the University of Edinburgh and earned me the *IMechE Best Project Award 2018*, a distinguished honor awarded once per university per year.


## Features
- **Social Force Model:** Implementation based on the model proposed by Helbing and Molnar.
- **Boarding Simulation:** Simulates the boarding process for an Airbus A320, generating boarding times and a newly-defined parameter, Time-To-Seat (TTS).
- **Comparison of Boarding Strategies:** Provides insights into various boarding strategies.

## My Contributions
- **MATLAB Code Development:** Created the MATLAB code and algorithm for modelling agents and simulating the boarding process.
- **New Metrics:** Developed the Time-To-Seat (TTS) parameter for comparing different boarding strategies.
- **Comparison of Boarding Strategies:** Provides insights into various boarding strategies, including a new method called RP3 (Reversed Pyramid 3)

The project includes simulations of various boarding methods. You can watch the simulations for each method at the following links:

- Wilma Method: https://youtu.be/DpsftdWtnVQ
- Random Method: https://youtu.be/gaFlylw3Hj0
- Block Method: https://youtu.be/jfdJho9hRuM
- RP3 Method: https://youtu.be/jM6nAA711D4

> [!NOTE]
> RP3 is over 20% better than boarding in blocks, where passengers board by sections from the back to the front (e.g., rows 30-40, then 20-29, and so on)

## Known Issues
- **Extensive Hard-Coding:** The code contains numerous hard-coded values, making it less flexible and harder to modify. Future iterations could benefit from parameterizing these values for greater adaptability and implementing OOP principles.
  ```
  if sex(i) == 1 && age(i) == 2
    ...
  elseif sex(i) == 1 && age(i) == 3
    ...
  elseif sex(i) == 1 && age(i) == 4
    ... 
  ```
    
- **Improvement Opportunities:** The flow and structure of the code can be optimized for better performance and maintainability.
