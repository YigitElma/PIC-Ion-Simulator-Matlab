classdef particle
    properties
        pos
        vel
        charge
        mass = 1
        trajectory
        track
    end

    methods
        function obj = particle(pos,vel,charge,mass,track)
            if nargin == 5
                obj.pos = pos;
                obj.vel = vel;
                obj.charge = charge;
                obj.mass = mass;
                obj.track = track;
            end
        end

        function obj = move(obj,F, dt)
            obj.vel = obj.vel + F/obj.mass * dt;
            obj.pos = obj.pos + obj.vel*dt; 

            if obj.track
                obj.trajectory = [obj.trajectory obj.pos'];
            end
        end

        function obj = MoveInField(obj,grad, dt)
            F = grad*obj.charge;
            obj.vel = obj.vel + F/obj.mass * dt;
            obj.pos = obj.pos + obj.vel*dt; 

            if obj.track
                obj.trajectory = [obj.trajectory obj.pos'];
            end
        end

        function obj = ReflectBottom(obj)
            obj.vel(2) = -obj.vel(2);
            obj.pos(2) = -obj.pos(2);
        end

        function obj = ReflectTop(obj,Lr)
            obj.vel(2) = -obj.vel(2);
            obj.pos(2) = Lr - (obj.pos(2) - Lr);
        end

        function obj = ReflectScreen(obj,x_screen)
            obj.pos(1) = x_screen - (obj.pos(1) - x_screen);
            obj.vel(1) = -obj.vel(1);
        end

        function obj = ReflectAccel(obj,x_accel)
            obj.pos(1) = x_accel - (obj.pos(1) - x_accel);
            obj.vel(1) = -obj.vel(1);
        end

    end
end
