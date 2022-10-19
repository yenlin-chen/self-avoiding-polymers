import numpy as np
from random import choice

class RandomWalkPolymer_xy():
    '''Builder for self-avoiding polymer chains using random walk.

    Object used for building self-avoiding polymer chains from 2D random
    walk, where the steps are restricted to either the x- or the y-
    direction. Diagonal steps are not allowed.
    '''

    def __init__(self, step_size=1):
        '''Class constructor.

        Params
            step_size (scalar): the size of one step

        Returns
            nothing
        '''
        self.step_size = step_size
        self.possible_steps = np.array([[-1,  0],
                                        [ 1,  0],
                                        [ 0, -1],
                                        [ 0,  1]]) * step_size

    def walk(self, current_pos):
        '''Returns a randomly selected position beside the current
        position.

        Params
            current_pos (ndarray): the current position

        Returns
            a position next to the input position base on random walk in
            the x- or y- axis

        '''
        # step randomly from current_pos position
        new_step = choice(self.possible_steps)
        return current_pos + new_step

    def get_neighbors(self, current_pos):
        '''Returns a list of all neighboring positions (4 positions).

        Params
            current_pos (ndarray): a position

        Returns
            list of neighbor positions
        '''
        neighbors_pos = []
        for move in self.possible_steps:
            neighbors_pos.append(current_pos+move)
        return neighbors_pos

    def is_restricted(self, pos_list):
        '''Checks if the current conformation is growable.

        Params
            pos_list (ndarray): list of positions

        Returns
            boolean value indicating whether conformation has reach its
            limit (True if yes)

        '''
        current_pos = pos_list[-1]

        neighbors_pos = self.get_neighbors(current_pos)

        cnt = 0
        for n in neighbors_pos:
            for p in pos_list:
                if np.array_equal(n, p):
                    cnt += 1
                    if cnt == 4:
                        return True
                    break

        return False

    def is_admissible(self, pos_list, next_pos):
        '''Checks if the new position is available for the given
        conformation.

        Params
            pos_list (ndarray): list of position for the current
                conformation
            next_pos (ndarray): candidate position for the next step

        Returns
            bool indicating whether the candidate position is valid for
            the given conformation
        '''
        return not any(np.array_equal(next_pos, p) for p in pos_list)

    def is_valid(self, pos_list):
        current_length = pos_list.shape[0]
        # nothing to check if there is only a starting position
        if current_length == 1:
            return True

        for idx in range(1, current_length):
            # check if length between node is equal to step size
            node_dist = np.linalg.norm(pos_list[idx]-pos_list[idx-1])
            if not node_dist == self.step_size:
                return False
            # check if shape is admissible
            if not self.is_admissible(pos_list[:idx], pos_list[idx]):
                return False

    def grow_to(self, target_length, pos_list=None, return_failed=False):
        '''Random walk until the chain reaches a specified length.

        Param
            target_length (int): length to reach
            pos_list (ndarray): list of position for the current
                conformation, the chain will grow from the last position
            return_failed (bool): whether to return a chain if its
                growth is restricted before reaching the specified
                target length

        Returns
            a chain of the target length grown by 2D random walk
        '''

        # start from origin if pos_list is None
        if pos_list is None:
            pos_list = np.zeros((target_length,2))
            length_to_grow = target_length - 1
        # else grow on pos_list
        else:
            # check if given conformation is valid
            if not self.is_valid(pos_list):
                raise ValueError('the given list of conformation is not valid')
            length_to_grow = target_length - pos_list.size
            pos_list = np.append(np.empty(length_to_grow,2))

        # iteratively generate positions
        for current_idx in range(length_to_grow):
            current_pos = pos_list[current_idx]

            if self.is_restricted(pos_list[:current_idx+1]):
                # print('fucked')
                if return_failed:
                    return pos_list[:current_idx+1]
                else:
                    try:
                        return self.grow_to(target_length)
                    except RecursionError as err:
                        print(err)
                        return False

            next_pos = self.walk(current_pos)
            while not self.is_admissible(pos_list[:current_idx+1], next_pos):
                next_pos = self.walk(current_pos)
            pos_list[current_idx+1] = next_pos

        return pos_list
