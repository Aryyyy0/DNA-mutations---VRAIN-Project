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
            score = inst.sample.log_acceptance_score
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