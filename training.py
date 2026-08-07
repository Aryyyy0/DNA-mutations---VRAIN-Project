"""
training.py - DNA Sequence Classification using 1D CNNs on Local Entropy Vectors
================================================================================
This module defines the PyTorch dataset structure for DNA sequence classification.
It processes in-memory instances of 'entropy_vector' containing local Shannon 
entropies and uses a log-likelihood threshold to partition typical vs. mutant samples.

Input Topology:
    - 1-D Grid of 10 sequential spatial coordinates (sequence windows)
    - Shape per sample: (1, 10) representing (in_channels, sequence_length)
"""

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
import probabilistic_automaton as pa
import torch.nn as nn

class DNAEntropyDataset(Dataset):
    def __init__(self, entropy_instances, threshold=-150.0):
        """
        Args:
            entropy_instances (list): List of 'entropy_vector' objects.
            threshold (float): Logarithmic probability threshold (theta).
                               Scores >= threshold are labeled 1 (Typical/Healthy).
                               Scores < threshold are labeled 0 (Low-Likelihood/Mutant).
        """
        self.vectors = []
        self.labels = []
        
        for inst in entropy_instances:
            # Extract the 10-dimensional local entropy vector
            self.vectors.append(inst.entropies)
            
            # Labelization: Evaluate the log-acceptance likelihood score against theta
            c = inst.sample.automata.accepts(inst.sample.sequence)
            score = c[1]
            label = 1 if score >= threshold else 0
            self.labels.append(label)           
        # Convert arrays to PyTorch-compatible float32 and int64 formats
        self.vectors = np.array(self.vectors, dtype=np.float32)
        self.labels = np.array(self.labels, dtype=np.int64)

    def __len__(self):
        return len(self.vectors)

    def __getitem__(self, idx):
        # Reshape the vector to (1, 10) representing (in_channels, sequence_length)
        # formatting is strictly expected by nn.Conv1d for 1-D convolutions.
        x = torch.tensor(self.vectors[idx]).unsqueeze(0)
        y = torch.tensor(self.labels[idx])
        return x, y
    
"""
Simple CNN for DNA Entropy Classification
 1-Stage convolutional block
 
"""
  
class DNAEntropyCNN_Simple(nn.Module):
    def __init__(self, num_classes=2):
        super(DNAEntropyCNN_Simple, self).__init__()
        
        # Calling  Conv1d with in_channels=1 because the input is a single-channel 1D vector (the entropy vector)
        # the output channels are set to 8, meaning we will learn 8 different filters
        self.conv = nn.Conv1d(in_channels=1, out_channels=8, kernel_size=3, padding=1)
        self.relu = nn.ReLU()
        self.pool = nn.AdaptiveAvgPool1d(1) # length reduction to 1, effectively summarizing the features across the sequence length
        
        # Output projection layer to categorical logits
        self.classifier = nn.Linear(in_features=8, out_features=2)

    def forward(self, x):
        # x has the format (Batch, 1, Longueur_sequence)
        
        x = self.conv(x)  
        x = self.relu(x)  
        x = self.pool(x)  
        
        x = x.squeeze(-1) # adapting the shape from (Batch, 8, 1) to (Batch, 8) for the linear layer
        
        logits = self.classifier(x) # -> (Batch, 2)
        return logits


"""
CNN for DNA Entropy Classification
 3-Stage convolutional blocks: Conv1d -> Non-linearity -> Pooling

coming soon ...

""" 