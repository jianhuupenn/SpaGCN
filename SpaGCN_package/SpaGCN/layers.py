import math

import torch
from torch.nn.modules.module import Module
from torch.nn.parameter import Parameter


class GraphConvolution(Module):
    """
    Simple GCN layer, similar to https://arxiv.org/abs/1609.02907
    """

    def __init__(
        self, in_features, out_features, bias=True, dtype=torch.float32, device="cpu"
    ):
        super(GraphConvolution, self).__init__()
        self.device = device
        self.dtype = dtype
        self.in_features = in_features
        self.out_features = out_features
        self.weight = Parameter(
            torch.rand(in_features, out_features, dtype=self.dtype, device=self.device)
        )
        if bias:
            self.bias = Parameter(
                torch.rand(out_features, dtype=self.dtype, device=self.device)
            )
        else:
            self.register_parameter("bias", None)
        self.reset_parameters()

    def reset_parameters(self):
        stdv = 1.0 / math.sqrt(self.weight.size(1))
        self.weight.data.uniform_(-stdv, stdv)
        if self.bias is not None:
            self.bias.data.uniform_(-stdv, stdv)

    def forward(self, input, adj):
        support = torch.mm(input, self.weight)
        output = torch.spmm(adj, support)
        if self.bias is not None:
            return output + self.bias
        else:
            return output

    def __repr__(self):
        return (
            self.__class__.__name__
            + " ("
            + str(self.in_features)
            + " -> "
            + str(self.out_features)
            + ")"
        )
