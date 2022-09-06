using System.Numerics;

namespace DelaunatorSharp
{
    public readonly struct Triangle
    {
        public readonly int Index;

        public readonly Vector2 Point0;

        public readonly Vector2 Point1;

        public readonly Vector2 Point2;


        public Triangle(int index, Vector2 point0, Vector2 point1, Vector2 point2)
        {
            Index = index;
            Point0 = point0;
            Point1 = point1;
            Point2 = point2;
        }
    }
}
